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


module iterative

  use ln_precision
  use ln_constants, only : pi
  use ln_allocation
  use mat_def
  use sparsekit_drv
  use inversions
  use ln_structure, only : TStruct_Info
  use lib_param, only : MAXNCONT, Tnegf, intarray, TInteractionArray
  use mpi_globals, only : id, numprocs, id0
  use outmatrix, only : outmat_c, inmat_c, direct_out_c, direct_in_c
  use clock
  !use transform

  implicit none
  private

  public :: calculate_transmissions
  public :: calculate_transmissions_and_dos

  public :: calculate_Gr
  public :: calculate_Gn_neq_components

  public :: rebuild_dns
  public :: calculate_gsmr_blocks
  public :: calculate_gsml_blocks
  public :: calculate_Gr_tridiag_blocks
  public :: calculate_Gr_column_blocks
  public :: calculate_Gn_tridiag_blocks
  public :: calculate_Gr_outer
  public :: calculate_Gn_outer

  public :: iterative_meir_wingreen
  public :: transmission_BP_corrected

  logical, parameter :: debug=.false.

  ! These are here temporarily ...
  type(z_DNS), dimension(:), allocatable :: gsmr
  type(z_DNS), dimension(:), allocatable :: gsml
  type(z_DNS), dimension(:,:), allocatable :: Gr
  type(z_DNS), dimension(:,:), allocatable :: ESH
  type(z_DNS), dimension(:,:), allocatable :: Gn


CONTAINS

  !****************************************************************************
  !
  ! Driver for computing Equilibrium Green Retarded (Outer blocks included)
  ! writing on memory
  !
  !****************************************************************************

  subroutine calculate_Gr(negf,E,SelfEneR,Tlc,Tcl,gsurfR,Grout,outer)

    !****************************************************************************
    !
    !Input
    !negf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR:   matrices array containing contacts surface green
    !outer:    optional parameter (0,1,2).
    !
    !Output:
    !Grout: Retarded Green's function (Device + Contacts overlap regions -> effective conductor)
    !   outer = 0  no outer parts are computed
    !   outer = 1  only D/C part is computed
    !   outer = 2  D/C and C/D parts are computed
    !              (needed for K-points calculations)
    !
    !*****************************************************************************

    type(Tnegf), intent(inout) :: negf
    complex(dp), intent(in) :: E
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    type(z_DNS), dimension(:), intent(inout) :: Tlc, Tcl, gsurfR
    type(z_CSR), intent(out) :: Grout
    integer, intent(in) :: outer

    !Work
    type(z_CSR) :: ESH_tot, Ain
    integer :: i,ierr, nbl, ncont,ii,n

    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts

    ! Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(negf%H,negf%S,(-1.0_dp, 0.0_dp),E,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot, ESH, negf%str%mat_PL_start)

    call destroy(ESH_tot)

    associate(cblk=>negf%str%cblk)
    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do
    end associate

    !! Add interaction self energy contribution, if any
    call add_sigma_r(negf%interactArray, ESH)

    call allocate_gsm(gsmr,nbl)
    call calculate_gsmr_blocks(ESH,nbl,2)

    call allocate_blk_dns(Gr,nbl)

    call calculate_Gr_tridiag_blocks(ESH,1)
    call calculate_Gr_tridiag_blocks(ESH,2,nbl)

    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    call set_Gr(negf%interactArray, Gr, negf%iE)

    call blk2csr(Gr,negf%str,negf%S,Grout)

    SELECT CASE (outer)
    CASE(0)
    CASE(1)
      call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.FALSE.,Grout)
    CASE(2)
      call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.TRUE.,Grout)
    end SELECT

    call destroy_blk(Gr)
    DEALLOCATE(Gr)

  end subroutine calculate_Gr

  !****************************************************************************
  !
  ! Driver for computing G_n contributions due to all contacts MINUS reference:
  ! Reference is necessary when splitting into contour + real-axis integration
  !
  !   Sum   [f_j(E)-f_r(E)] Gr Gam_j Ga
  !   j!=r
  !
  ! NOTE: Setting reference to ncont+1 removes the reference part
  !****************************************************************************

  subroutine calculate_Gn_neq_components(negf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,Glout,outblocks)

    !****************************************************************************
    !
    !Input
    !negf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR:   matrices array containing contacts surface green
    !outblocks: optional parameter (0,1,2).
    !
    !Output:
    !Gn: NE GF (Device + Contacts overlap regions -> effective conductor)
    !   outblocks = 0  no outer parts are computed
    !   outblocks = 1  only D/C part is computed
    !   outblocks = 2  D/C and C/D parts are computed
    !              (needed for K-points calculations)
    !
    !*****************************************************************************

    implicit none

    !In/Out
    type(Tnegf), intent(inout) :: negf
    type(z_CSR), intent(inout)  :: Glout
    type(z_DNS), dimension(:), intent(in)  :: SelfEneR, gsurfR, Tlc, Tcl
    real(dp), intent(in)  :: E
    real(dp), dimension(:), intent(in)  :: frm
    integer, intent(in)  :: outblocks

    !Work
    integer :: ref
    complex(dp) :: Ec
    integer :: i,ierr,ncont,nbl, lbl, rbl
    integer, dimension(:), allocatable :: Gr_columns
    type(z_CSR) :: ESH_tot, Gl
    logical :: mask(MAXNCONT)


    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts
    ref = negf%refcont
    associate (cblk=>negf%str%cblk, indblk=>negf%str%mat_PL_start)

    Ec = cmplx(E,0.0_dp,dp)

    ! Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(negf%H,negf%S,(-1.0_dp, 0.0_dp),Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot,ESH, negf%str%mat_PL_start)

    call destroy(ESH_tot)

    ! Add contact self-energy to proper blocks
    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do

    ! Add interaction self energy if any and initialize scba counter
    call add_sigma_r(negf%interactArray, ESH)

    call allocate_gsm(gsmr,nbl)
    call allocate_gsm(gsml,nbl)

    ! Determine the leftmost and rightmost contact blocks to determine
    ! which column blocks are needed and hence which gsmr and gsml. In the case
    ! of interactions we need to be able to calculate all columns.
    if (allocated(negf%interactArray)) then
      rbl = 1
      lbl = nbl - 1
    else
      mask = .true.
      mask(ref) = .false.
      rbl = minval(cblk(1:ncont),mask(1:ncont))
      lbl = maxval(cblk(1:ncont),mask(1:ncont))
    endif

    call calculate_gsmr_blocks(ESH,nbl,rbl+1)
    call calculate_gsml_blocks(ESH,1,lbl-1)

    call allocate_blk_dns(Gr,nbl)

    ! -------------------------------------------------------------
    ! 1. rbl>lbl  => lbl+1=rbl-1 => compute first Gr(rbl-1,rbl-1)
    ! 2. rbl<lbl  => lbl=rbl-2 has been computed
    ! calculate_Gr does not compute if sbl>nbl or sbl<1
    call calculate_Gr_tridiag_blocks(ESH,rbl)
    call calculate_Gr_tridiag_blocks(ESH,rbl+1,nbl)
    call calculate_Gr_tridiag_blocks(ESH,rbl-1,1)
    !Passing Gr to interaction that builds Sigma_n
    call set_Gr(negf%interactArray, Gr, negf%iE)

    !Computes the columns of Gr for the contacts != reference
    ! Keep track of the calculated column indices in the array Gr_columns.
    allocate(Gr_columns(ncont))
    Gr_columns = 0
    do i=1,ncont
      if (i.NE.ref) THEN
        Gr_columns(i) = cblk(i)
        call calculate_Gr_column_blocks(ESH,cblk(i),indblk)
      endif
    end do

    !If interactions are not present we can already destroy gsmr, gsml.
    !Otherwise they are still needed to calculate columns on-the-fly.
    if (.not.allocated(negf%interactArray)) then
      call destroy_gsm(gsmr)
      call deallocate_gsm(gsmr)
      call destroy_gsm(gsml)
      call deallocate_gsm(gsml)
    end if

    ! Computing device G_n
    call allocate_blk_dns(Gn,nbl)
    call init_blkmat(Gn,ESH)
    ! First the coherent Gn is computed from contact self-energies
    call calculate_Gn_tridiag_blocks(ESH,SelfEneR,frm,ref,negf%str,Gn)

    ! Adding el-ph part: G^n = G^n + G^r Sigma^n G^a (at first call does nothing)
    ! NOTE:  calculate_Gn has factor [f_i - f_ref], hence all terms will contain this factor
    if (allocated(negf%interactArray)) then
      call calculate_Gn_tridiag_elph_contributions(negf,ESH,Gn,Gr_columns)
      call destroy_gsm(gsmr)
      call deallocate_gsm(gsmr)
      call destroy_gsm(gsml)
      call deallocate_gsm(gsml)
    end if


    !Passing G^n to interaction that builds Sigma^n
    call set_Gn(negf%interactArray, Gn, negf%iE)

    call blk2csr(Gn,negf%str,negf%S,Glout)

    end associate

    !Computing the 'outer' blocks (device/contact overlapping elements)
    SELECT CASE (outblocks)
    CASE(0)
    CASE(1)
      call calculate_Gn_outer(Tlc,gsurfR,SelfEneR,negf%str,frm,ref,.false.,Glout)
    CASE(2)
      call calculate_Gn_outer(Tlc,gsurfR,SelfEneR,negf%str,frm,ref,.true.,Glout)
    end SELECT

    if (negf%tDestroyBlk) then
      call destroy_all_blk()
    end if

  end subroutine calculate_Gn_neq_components

  !---------------------------------------------------------------------
  !>
  !  Iterative algorithm implementing Meir Wingreen formula for a given
  !  electrode
  !  Note: self consistent born approximation is not accounted for here
  !  It is assumed that the el-ph container already includes the
  !  desired values. SCBA loop should be run outside
  !  Many operations from calls_neq_ph are repeated here, as it is
  !  assumed that A and Gn are not available at the time of the call
  !
  !  It implements the form without the reference electrode:
  !
  !  I_i = Tr[\Sigma_{i}^{n}A - \Gamma_{i}G^{n}] =
  !      = Tr[\Gamma_{i}(f_{i}A - f_{ref}A - G^{n,l\neq ref} - G^{n,l\neq ref}_{\phi})]
  !
  !  where G^{n,l\neq ref} is the component including no el-ph
  !  If i=ref it reduces to
  !  I_i = Tr[ \Gamma_{i}(-G^{n,l\neq ref}  - G^{n,l\neq ref}_{\phi}) ]
  !---------------------------------------------------------------------
  subroutine iterative_meir_wingreen(negf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,curr_mat)
    type(Tnegf), intent(inout) :: negf
    real(dp) :: E
    type(z_DNS), dimension(:) :: SelfEneR, gsurfR, Tlc, Tcl
    real(dp), dimension(:) :: frm
    real(dp), dimension(:) :: curr_mat

    !Work
    complex(dp) :: Ec, tmp
    integer :: i,ierr,ncont,nbl,lbl,ii,n
    integer :: pl_start, pl_end
    integer :: ref, lead, lead_blk, ref_blk
    integer, dimension(:), allocatable :: Gr_columns
    type(z_DNS) :: work1, work2, Gam, A
    type(z_CSR) :: ESH_tot, Gl

    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts
    ref = negf%refcont
    ref_blk = negf%str%cblk(ref)
    associate (cblk=>negf%str%cblk)

    Ec=cmplx(E,0.0_dp,dp)

    call prealloc_sum(negf%H,negf%S,(-1.0_dp, 0.0_dp),Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot, ESH, negf%str%mat_PL_start)

    call destroy(ESH_tot)

    ! Add contact self energies
    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do

    !! Add el-ph self energy if any
    call add_sigma_r(negf%interactArray, ESH)

    call allocate_gsm(gsmr,nbl)
    call allocate_gsm(gsml,nbl)

    call calculate_gsmr_blocks(ESH,nbl,2)

    if ( ncont.gt.1 ) then
      call calculate_gsml_blocks(ESH,1,nbl-2)
    endif

    call allocate_blk_dns(Gr,nbl)

    call calculate_Gr_tridiag_blocks(ESH,1)
    call calculate_Gr_tridiag_blocks(ESH,2,nbl)

    !! Give Gr to interaction model if any
    call set_Gr(negf%interactArray, Gr, negf%iE)

    !! Never calculate outer blocks
    call allocate_blk_dns(Gn, nbl)

    call init_blkmat(Gn,ESH)

    !! TEMPORARY AND INEFFICIENT:
    !! CALCULATE THE FULL Gn WHEN THE CONTACT BLOCK WHOULD BE ENOUGH
    allocate(Gr_columns(ncont))
    Gr_columns = 0
    do i=1,ncont
      if (i.NE.ref) THEN
        Gr_columns(i) = cblk(i)
        call calculate_Gr_column_blocks(ESH, cblk(i), negf%str%mat_PL_start)
      endif
    end do
    call calculate_Gn_tridiag_blocks(ESH, SelfEneR, frm, ref, negf%str, Gn)

    if (allocated(negf%interactArray)) then
      call calculate_Gn_tridiag_elph_contributions(negf, ESH, Gn, Gr_columns)
    end if

    ! The gsmr, gsml are used to calculate columns on-the-fly and now can be removed
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)
    call destroy_gsm(gsml)
    call deallocate_gsm(gsml)

    do i=1,size(negf%ni)
      lead = negf%ni(i)
      lead_blk = cblk(lead)
      call zspectral(SelfEneR(lead),SelfEneR(lead), 0, Gam)
      if (lead.eq.ref) then
        call prealloc_mult(Gam, Gn(lead_blk, lead_blk), (-1.0_dp, 0.0_dp), work1)
        curr_mat(i) = real(trace(work1))
        call destroy(work1)
      else
        call zspectral(Gr(lead_blk, lead_blk), Gr(lead_blk, lead_blk), 0, A)
        tmp = frm(lead)-frm(ref)
        call prealloc_sum(A, Gn(lead_blk, lead_blk), tmp, (-1.0_dp, 0.0_dp), work1)
        call destroy(A)
        call prealloc_mult(Gam, work1, work2)
        call destroy(work1)
        curr_mat(i) = real(trace(work2))
        call destroy(work2)
      endif
      call destroy(Gam)
    end do
    !Convert to output CSR format.
    call blk2csr(Gn,negf%str,negf%S,Gl)

    call destroy_all_blk()

    end associate

  end subroutine iterative_meir_wingreen

  !---------------------------------------------------------------------
  !  Calculate the layer current per unit energy
  !
  !    I_LL'(E) = Tr[(ES-H)_LL' * Gn_L'L(E)-(ES-H)_L'L * Gn_LL'(E)]
  !
  subroutine iterative_layer_current(negf,E,curr_mat)
    type(Tnegf), intent(inout) :: negf
    real(dp) :: E
    real(dp), dimension(:) :: curr_mat

    integer :: nbl, ii
    type(z_DNS) :: work1
    complex(dp), parameter :: minusone = (-1.0_dp, 0.0_dp)

    nbl = negf%str%num_PLs

    if (size(curr_mat) .ne. nbl-1) then
         stop 'ERROR: curr_mat with wrong size in iterative_layer_current'
   end if

    do ii = 1, nbl-1
      call prealloc_mult(ESH(ii,ii+1),Gn(ii+1,ii),work1)
      call prealloc_mult(ESH(ii+1,ii),Gn(ii,ii+1), minusone, work1)
      curr_mat(ii) = trace(work1) 
    end do
    
    if (negf%tDestroyBlk) then
      call destroy_all_blk()
    end if

  end subroutine iterative_layer_current


  !---------------------------------------------------------------------
  subroutine destroy_all_blk()
    call destroy_blk(Gr)
    DEALLOCATE(Gr)
    call destroy_blk(Gn)
    DEALLOCATE(Gn)
    call destroy_ESH(ESH)
    DEALLOCATE(ESH)
  end subroutine

  ! Add all Self-enrgies to the Hamiltonian
  subroutine add_sigma_r(interactArray, ESH)
    type(TInteractionArray), allocatable :: interactArray(:)
    type(z_DNS) :: ESH(:,:)

    integer :: kk
    if (allocated(interactArray)) then
      do kk = 1, size(interactArray)
        call interactArray(kk)%inter%add_sigma_r(ESH)
      end do
    end if
  end subroutine add_sigma_r

  ! Provides Gr to the interaction models.
  ! In some case the self/energies are computed
  subroutine set_Gr(interactArray, Gr, iE)
    type(TInteractionArray), allocatable :: interactArray(:)
    type(z_DNS) :: Gr(:,:)
    integer :: iE

    integer :: kk
    if (allocated(interactArray)) then
      do kk = 1, size(interactArray)
        call interactArray(kk)%inter%set_Gr(Gr, iE)
      end do
    end if
  end subroutine set_Gr

  ! Provides Gn to the interaction models.
  ! In some case the self/energies are computed
  subroutine set_Gn(interactArray, Gn, iE)
    type(TInteractionArray), allocatable :: interactArray(:)
    type(z_DNS) :: Gn(:,:)
    integer :: iE

    integer :: kk
    if (allocated(interactArray)) then
      do kk = 1, size(interactArray)
        call interactArray(kk)%inter%set_Gn(Gn, iE)
      end do
    end if
  end subroutine set_Gn

  !------------------------------------------------------------------------------!
  ! Transmission_BP_corrected
  !
  ! The routine implements the current-probes, i.e. it assumes that the BPs are
  ! elastic dephasing sources that fullfill current conservation at any energy
  ! For a nice discussion see papers by D. Ryndyk and
  ! M. Kilgour, D. Segal, The J. of Chem Phys 144, 124107 (2016)
  !
  ! LIMITATIONS:
  !
  ! This routine is limited to 2 contacts
  ! It uses dense matrices and assumes there is only 1 PL in the device
  !------------------------------------------------------------------------------!

  subroutine transmission_BP_corrected(negf,SelfEneR,tun_mat)

    type(Tnegf), intent(inout) :: negf
    type(z_DNS), dimension(:) :: SelfEneR
    real(dp), dimension(:) :: tun_mat

    integer :: nn, mm, cont, ni, nf
    integer :: NumOrbs, ncont
    integer, allocatable ::  pls(:), ple(:)
    real(dp), allocatable  :: Trans(:,:), W(:,:), GG(:,:), R(:)
    type(z_DNS) :: Gam1, Gam2, GreenR, GreenA, Tmp1, Tmp2


    NumOrbs=negf%str%central_dim
    ncont = negf%str%num_conts

    allocate(Trans(NumOrbs+ncont,NumOrbs+ncont))
    Trans=0.0_dp

    call create(Gam1, NumOrbs, NumOrbs)
    call create(Gam2, NumOrbs, NumOrbs)

    call blk2dns(Gr,negf%str,GreenR)
    call zdagger(GreenR, GreenA)

    allocate(pls(ncont))
    allocate(ple(ncont))

    do nn = 1, ncont
      pls(nn) = negf%str%mat_PL_start(negf%str%cblk(nn))
      ple(nn) = negf%str%mat_PL_end(negf%str%cblk(nn))
    end do

    do nn = 1, NumOrbs+ncont

      if (nn > NumOrbs) then
         cont = nn - NumOrbs
         call zspectral(SelfEneR(cont),SelfEneR(cont),0,Tmp1)
         Gam1%val(pls(cont):ple(cont),pls(cont):ple(cont)) = Tmp1%val
         call destroy(Tmp1)
      else
         Gam1%val = 0.0_dp
         Gam1%val(nn,nn)=negf%bp_deph%coupling(nn)
      end if

      do mm = 1, NumOrbs+ncont

        if (mm > NumOrbs) then
           cont = mm - NumOrbs
           call zspectral(SelfEneR(cont),SelfEneR(cont),0,Tmp2)
           Gam2%val(pls(cont):ple(cont),pls(cont):ple(cont)) = Tmp2%val
           call destroy(Tmp2)
        else
           Gam2%val = 0.0_dp
           Gam2%val(mm,mm)=negf%bp_deph%coupling(mm)
        end if

        ! Compute coherent transmission: Tr[Gam1 Gr Gam2 Ga]
        ! The inefficient quadruple loop has been substituted with
        ! M*M multiplications exploting OMP parallelism
        call prealloc_mult(Gam1,GreenR,Tmp1)
        call prealloc_mult(Tmp1,Gam2,Tmp2)
        call destroy(Tmp1)
        call prealloc_mult(Tmp2,GreenA,Tmp1)
        call destroy(Tmp2)

        Trans(nn,mm) = real(trace(Tmp1))

        call destroy(Tmp1)

      end do
    end do

    call destroy(Gam1, Gam2, GreenR, GreenA)
    deallocate(pls, ple)

    allocate(W(NumOrbs,NumOrbs))
    allocate(R(NumOrbs+ncont))

    do nn = 1, NumOrbs+ncont
       R(nn) = 1.0_dp
       do mm= 1, NumOrbs+ncont
          if (mm.ne.nn) then
             R(nn) = R(nn)-Trans(mm,nn)
          end if
       end do
    end do

    do nn = 1, NumOrbs
      do mm = 1, NumOrbs
         W(nn,mm) = -Trans(nn,mm)
      end do
      W(nn,nn) = 1.0_dp - R(nn)
    end do

    deallocate(R)
    allocate(GG(NumOrbs,NumOrbs))

    call inverse(GG,W,NumOrbs)

    deallocate(W)

    allocate(R(NumOrbs))

    do nn = 1, size(negf%ni)
      ni = negf%ni(nn)
      nf = negf%nf(nn)
      R =  matmul(Trans(NumOrbs+ni,:),GG)
      tun_mat(nn) = Trans(NumOrbs+ni,NumOrbs+nf) + dot_product(R, Trans(:,NumOrbs+nf))
    end do

    deallocate(Trans,GG)
    deallocate(R)

  end subroutine transmission_BP_corrected

  !------------------------------------------------------------------------------!
  !DAR end
  !------------------------------------------------------------------------------!

  !**********************************************************************
  subroutine init_blkmat(Matrix,S)
    type(z_DNS), dimension(:,:) :: Matrix,S

    integer :: nbl, j

    nbl = SIZE(Matrix,1)

    call create(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)
    Matrix(1,1)%val=(0.0_dp,0.0_dp)
    do j=2,nbl-1
      call create(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)
      Matrix(j-1,j)%val=(0.0_dp,0.0_dp)
      call create(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)
      Matrix(j,j)%val=(0.0_dp,0.0_dp)
      call create(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)
      Matrix(j,j-1)%val=(0.0_dp,0.0_dp)
    end do
    if (nbl.gt.1) then
      call create(Matrix(nbl,nbl),S(nbl,nbl)%nrow,S(nbl,nbl)%ncol)
      Matrix(nbl,nbl)%val=(0.0_dp,0.0_dp)
      call create(Matrix(nbl-1,nbl),S(nbl-1,nbl)%nrow,S(nbl-1,nbl)%ncol)
      Matrix(nbl-1,nbl)%val=(0.0_dp,0.0_dp)
      call create(Matrix(nbl,nbl-1),S(nbl,nbl-1)%nrow,S(nbl,nbl-1)%ncol)
      Matrix(nbl,nbl-1)%val=(0.0_dp,0.0_dp)
    endif

  end subroutine init_blkmat


  !**********************************************************************
  subroutine destroy_gsm(gsm)
    type(z_DNS), dimension(:) :: gsm
    integer :: i, i1, nbl

    nbl=size(gsm,1)

    do i=1,nbl
      if (allocated(gsm(i)%val)) call destroy(gsm(i))
    end do

  end subroutine destroy_gsm

  !**********************************************************************
  subroutine destroy_blk(M)
    type(z_DNS), dimension(:,:) :: M
    integer :: i, i1, nbl

    nbl=size(M,1)

    do i=1,nbl
      do i1=1,nbl
        if (ALLOCATED(M(i1,i)%val)) THEN
          call destroy(M(i1,i))
        end if
      end do
    end do

  end subroutine destroy_blk

  !**********************************************************************
  subroutine destroy_ESH(ESH)

    integer :: i, nbl
    type(z_DNS), dimension(:,:) :: ESH

    nbl=size(ESH,1)

    do i=1,nbl
      call destroy(ESH(i,i))
    end do
    do i=2,nbl
      call destroy(ESH(i-1,i))
      call destroy(ESH(i,i-1))
    end do

  end subroutine destroy_ESH


  !***********************************************************************
  !
  !  Construct a sparse matrix starting from the sparse matrices array
  !  related to the blocks
  !
  !***********************************************************************

  subroutine rebuild_dns(Atot,A,n,indb)

    !***********************************************************************
    !Input:
    !A: sparse matrices array
    !n: A dimension (A(n,n))
    !indb: blocks indeces array (contains position of first elements of diagonal
    !      blocks, included an additional ending address "last row+1")
    !
    !Output:
    !Atot: total sparse matrix (allocated internally)
    !***********************************************************************

    implicit none

    !In/Out
    type(z_CSR) :: Atot  !Allocato internamente
    integer :: n
    integer, dimension(n+1) :: indb
    type(z_cSR), dimension(n,n) :: A

    !Work
    integer :: i,j,Atot_nrow,i1,j1

    Atot_nrow=indb(n+1)-1

    call create(Atot,Atot_nrow,Atot_nrow,0)
    Atot%rowpnt(:)=1

    do i=1,n
      do j=1,n

        if (A(i,j)%nrow.GT.0) THEN
          i1=indb(i)
          j1=indb(j)
          call concat(Atot,A(i,j),i1,j1)
        endif

      end do
    end do

  end subroutine rebuild_dns

  !***********************************************************************
  !
  !  g_small right (gsmr) calculation - write on memory
  !
  !***********************************************************************

  subroutine calculate_gsmr_blocks(ESH,sbl,ebl,keepall)

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




  !***********************************************************************
  !
  !  g_small left (gsml) calculation - write on memory
  !
  !***********************************************************************

  subroutine calculate_gsml_blocks(ESH,sbl,ebl)

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

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'calculate_gsml done'
      WRITE(*,*) '********************'
    endif

  end subroutine calculate_gsml_blocks





  !***********************************************************************
  !
  !  Diagonal, Subdiagonal, Superdiagonal blocks of Green Retarded
  !  Gr(nbl,nbl) - writing on memory
  !
  !***********************************************************************

  subroutine calculate_Gr_tridiag_blocks(ESH,sbl,ebl)

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
    type(z_DNS), dimension(:,:) :: ESH
    integer :: sbl
    integer, optional :: ebl

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

  !**************************************************************************
  !
  !  Calculate Green Retarded column "n" - writing on memory
  !
  !**************************************************************************

  subroutine calculate_Gr_column_blocks(ESH,n,indblk)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !n: n umber of column to be calculated
    !
    !global variables needed: nbl (number of layers), indblk(nbl+1),
    !Gr diagonal, subadiagonal and superdiagonal, gsmr(:) for
    !downgoing and gsml(:) for upgoing
    !
    !Output:
    !sparse matrices array global variable Gr(:,n) is available in
    !memory - single blocks are allocated internally, array Gr(nbl,nbl)
    !must be allocated externally
    !***********************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: n
    integer, dimension(:), intent(in) :: indblk

    !Work
    integer :: i,nrow,ncol,nbl
    type(z_DNS) :: work1
    real(dp) :: max

    nbl = size(ESH,1)

    if (n.GT.nbl) THEN
      STOP 'Error in calculate_Grcol : n is greater than nbl'
    endif

    !***************************************
    !  Downgoing (j>=n+2 && n<nbl-1)
    !
    !   G_j,n = -gR_jj T_j,j-1 G_j-1,n
    !
    !***************************************
    if (n.LT.(nbl-1)) THEN

      do i=n+2,nbl

        max=MAXVAL(ABS(Gr(i-1,n)%val))
        if (max.GT.EPS) THEN
          call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.0_dp, 0.0_dp),work1)
          call prealloc_mult(work1,Gr(i-1,n),Gr(i,n))
          call destroy(work1)
        else
          ! WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
          exit
        end if

      end do

    endif
    !*************************************
    !   Up-going (j<=n-2 && n>2)
    !
    !   G_j,n = -gL_jj T_j,j+1 G_j+1,n
    !
    !*************************************

    if (n.GT.2) THEN

      do i=n-2,1,(-1)

        max=MAXVAL(ABS(Gr(i+1,n)%val))

        if (max.GT.EPS) THEN
          call prealloc_mult(gsml(i),ESH(i,i+1),(-1.0_dp, 0.0_dp),work1)
          call prealloc_mult(work1,Gr(i+1,n),Gr(i,n))
          call destroy(work1)
        else
          ! WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
          exit
        endif

      end do

    endif

    if (debug) then
      WRITE(*,*) '******************************'
      WRITE(*,*) 'calculate_Grcol done column',n
      WRITE(*,*) '******************************'
    endif

  end subroutine calculate_Gr_column_blocks

  !****************************************************************************
  !
  !  Calculate Green Retarded - writing on memory (optimized on mask)
  !
  !****************************************************************************

  subroutine Gr_blk2csr(P,nbl,indblk,A)

    !****************************************************************************
    !Input:
    !P: CSR matrix containing masking pattern
    !
    !global variable needed: nbl, indblk(nbl+1), Gr(:,:)
    !
    !Output:
    !A: sparse matrix containing Green Retarded of device (allocated internally)
    !****************************************************************************

    implicit none

    !In/Out
    integer :: nbl
    integer, dimension(:), pointer :: indblk
    type(z_CSR) :: A, P, GrCsr

    !Work
    integer :: i, j, i1, ix, iy, x, y, col, oldx

    !create A with same pattern of P
    call create(A,P%nrow,P%ncol,P%nnz)
    A%rowpnt(:)=P%rowpnt(:)
    A%colind(:)=P%colind(:)
    A%nzval = (0.0_dp,0.0_dp)

    !If only one block is present, concatenation is not needed
    !and it's implemented in a more trivial way
    if (nbl.EQ.1) THEN

      call create(GrCsr,Gr(1,1)%nrow,Gr(1,1)%ncol,Gr(1,1)%nrow*Gr(1,1)%ncol)
      call dns2csr(Gr(1,1),GrCsr)

      call mask(GrCsr,P,A)
      call destroy(GrCsr)

    ELSE

      !Cycle upon all rows
      x = 1
      do i = 1, A%nrow
        !Choose which block (row) we're dealing with
        oldx = x

        !Check if row is in same block of previous or in next block. Not needed
        !(and not allowed not to exceed indblk index boundaries) if we're in the last block
        if (oldx.EQ.nbl) THEN
          x = oldx
        ELSE
          do ix = oldx, oldx+1
            if ( (i.GE.indblk(ix)).AND.(i.LT.indblk(ix+1)) ) x = ix
          end do
        endif

        !Offset: i1 is the index for separate blocks
        i1 = i - indblk(x) + 1
        !Cycle upon columns
        do j = A%rowpnt(i), A%rowpnt(i+1) -1
          !Choose which block column we're dealing with
          y = 0
          if (x.EQ.1) THEN
            if ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then
              y = 1
            ELSEif ( (A%colind(j).GE.indblk(x + 1)).AND.(A%colind(j).LT.indblk(x + 2)) ) then
              y = 2
            endif
          elseif (x.eq.nbl) then
            if ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then
              y = nbl
            ELSEif ( (A%colind(j).GE.indblk(x - 1)).AND.(A%colind(j).LT.indblk(x)) ) then
              y = nbl - 1
            endif
          ELSE
            do iy = x-1, x+1
              if ( (A%colind(j).GE.indblk(iy)).AND.(A%colind(j).LT.indblk(iy + 1)) ) y = iy
            end do
          endif
          if (y.eq.0) then
            write(*,*)
            write(*,*) 'ERROR in blk2csr: probably wrong PL size',x
            write(*,*) 'row',i,A%colind(j)
            write(*,*) 'block indeces:',indblk(1:nbl)
            stop
          endif

          col = A%colind(j) - indblk(y) + 1

          A%nzval(j) = Gr(x,y)%val(i1,col)

        end do

      end do

    endif

    !if (debug) call writePeakInfo(6)
    if (debug) then
      WRITE(*,*) '**********************'
      WRITE(*,*) 'calculate_GreenR done'
      WRITE(*,*) '**********************'
    endif

  end subroutine Gr_blk2csr


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

  subroutine calculate_Gn_tridiag_blocks(ESH,SelfEneR,frm,ref,struct,Gn)

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
    type(z_DNS), dimension(:,:), intent(in) :: ESH
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

  !****************************************************************************
  !
  ! Calculate G_n contributions for all contacts (except reference)
  ! This version computes Grcol on the fly
  !
  !****************************************************************************

  subroutine calculate_Gn_tridiag_blocks2(ESH,SelfEneR,frm,ref,struct,Gn)

    !******************************************************************************
    !Input:
    !ESH(nbl,nbl): sparse matrices array ES-H
    !SelfEneR(ncont): sparse matrices array with contacts Self Energy
    !frm(ncont): Fermi diistribution value for each contact
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
    type(z_DNS), dimension(:,:), intent(in) :: ESH
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
        if (Gr(1,cb)%nrow.gt.0) then
          call prealloc_mult(Gr(1,cb),Gam,work1)
          call zdagger(Gr(1,cb),Ga)
          call prealloc_mult(work1,Ga,frmdiff,Gn(1,1))
          call destroy(work1, Ga)
        else
          Gn(1,1)%val=(0.0_dp,0.0_dp)
        endif

        ! Computation of all tridiagonal blocks
        do i=2,nbl

          ! Computation of Gr(i,cb) assuming Gr(i-1,cb) exists
          ! Assume downgoing: i > cb
          if (Gr(i-1,cb)%nrow.GT.0) THEN
            call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.0_dp, 0.0_dp),work1)
            call destroy(gsmr(i))
            call prealloc_mult(work1,Gr(i-1,cb),Gr(i,cb))
            call destroy(work1)
            if (MAXVAL(ABS(Gr(i,cb)%val)).lt.EPS) call destroy(Gr(i,cb))
          endif

          ! Computation of Gl(i,j) = Gr(i,cb) Gam(cb) Ga(cb,j)
          ! Both Gr(i,cb) and Gr(j,cb) must be non-zero
          if (Gr(i-1,cb)%nrow.gt.0 .and. Gr(i,cb)%nrow.gt.0) THEN
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
          ELSE
            Gn(i  ,i)%val=(0.0_dp,0.0_dp)
            Gn(i-1,i)%val=(0.0_dp,0.0_dp)
            Gn(i,i-1)%val=(0.0_dp,0.0_dp)
          endif

          if (Gr(i-1,cb)%nrow.gt.0) call destroy(Gr(i-1,cb))

        end do

        call destroy(Gam)

      endif

    end do

  end subroutine calculate_Gn_tridiag_blocks2

  !****************************************************************************
  !
  ! Calculate G_n contributions due to elph:  G_n = G_n + Gr Sigma_ph Ga
  !
  !****************************************************************************
  subroutine calculate_Gn_tridiag_elph_contributions(negf, ESH, Gn, existing_Gr_cols)

    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    type(z_DNS), dimension(:,:), intent(inout) :: Gn
    integer, dimension(:), intent(in) :: existing_Gr_cols

    Type(z_DNS), dimension(:,:), allocatable :: Sigma_ph_n, sigma_blk
    Type(z_DNS) :: Ga, work1, work2, sigma_tmp
    integer :: n, k, nbl, nrow, ierr, ii, jj, norbs, nblk, indstart, indend

    nbl = negf%str%num_PLs

    !! If this is the first scba cycle, there's nothing to do
    if (negf%scbaDriver%scba_iter .eq. 0) then
      return
    endif
    allocate(Sigma_ph_n(nbl,nbl),stat=ierr)
    if (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'

    ! The block sigma n is made available from el-ph model
    ! Note: the elph models could not keep a copy and calculate it
    ! on the fly. You have to rely on the local copy
    if (allocated(negf%interactArray)) then
      call negf%interactArray(1)%inter%get_sigma_n(Sigma_ph_n, negf%ie)
    end if

    !! Calculate the diagonal and off diagonal (if needed) blocks of Gn
    !! in the assumption of diagonal self energy
    !! Gn(n,n) = Gr(n,k)*Sigma_n(k,k)*Ga(k,n)
    !! Gn(n,n+1) = Gr(n,k)*Sigma_n(k,k)*Ga(k,n+1)
    !! Gn(n,n-1) = Gr(n,k)*Sigma_n(k,k)*Ga(k,n-1)
    !! All the rows of Gr need to be available

    do k = 1, nbl
      ! Calculate k-th column on-the-fly. The reference contact column
      ! might already be available. Check if the top and bottom of the
      ! column are available.
      if (all(existing_Gr_cols .ne. k)) then
        call calculate_Gr_column_blocks(ESH, k, negf%str%mat_PL_start)
      endif

      do n = 1, nbl

        ! Zero column blocks are not created at all.
        if (allocated(Gr(n, k)%val)) then
          ! Calculate diagonal blocks Gn(n, n)
          call zdagger(Gr(n, k), Ga)
          call prealloc_mult(Gr(n, k), Sigma_ph_n(k, k), work1)
          call prealloc_mult(work1, Ga, work2)
          Gn(n,n)%val = Gn(n,n)%val + work2%val
          call destroy(work2, Ga)

          ! Computing blocks of Gn(n, n - 1).
          if (n .lt. nbl .and. allocated(Gn(n, n + 1)%val) .and. allocated(Gr(n + 1, k)%val)) then
            call zdagger(Gr(n + 1, k), Ga)
            call prealloc_mult(work1, Ga, work2)
            Gn(n, n + 1)%val = Gn(n, n + 1)%val + work2%val
            call destroy(work2, Ga)
          endif

          ! Computing blocks of Gn(n, n - 1).
          if (n .gt. 1 .and. allocated(Gn(n, n - 1)%val) .and. allocated(Gr(n - 1, k)%val)) then
            call zdagger(Gr(n - 1, k), Ga)
            call prealloc_mult(work1, Ga, work2)
            Gn(n, n - 1)%val = Gn(n, n - 1)%val + work2%val
            call destroy(work2, Ga)
          endif
          call destroy(work1)
        end if
      end do

      ! Remove column blocks of Gr if they were calculated on the fly.
      do n = 1, nbl
        if (abs(n - k) .gt. 1 .and. all(existing_Gr_cols .ne. k)) then
          call destroy(Gr(n, k))
        end if
      end do

    end do

    deallocate(Sigma_ph_n)

  end subroutine calculate_Gn_tridiag_elph_contributions

  subroutine blk2dns(G,str,Gdns)
    type(z_DNS), dimension(:,:), intent(in) :: G
    type(Tstruct_info), intent(in) :: str
    type(z_DNS) :: Gdns

    integer :: n, m, ii, jj, pl_start1, pl_end1, pl_start2, pl_end2, nbl, nrows

    nbl = str%num_PLs
    nrows = str%central_dim
    call create(Gdns,nrows,nrows)

    do n = 1, nbl
      do m= 1, nbl
         pl_start1 = str%mat_PL_start(n)
         pl_end1 = str%mat_PL_end(n)
         pl_start2 = str%mat_PL_start(m)
         pl_end2 = str%mat_PL_end(m)

         do ii = 1, pl_end1 - pl_start1 + 1
           do jj = 1, pl_end2 - pl_start2 + 1
               Gdns%val(pl_start1 + ii - 1,pl_start2 + jj - 1) = G(n,m)%val(ii,jj)
           end do
         end do
      end do
    end do

  end subroutine blk2dns


  !Concatenation for every contact in G_n. Performs a sum on elements, not a replacement
  !Similar to calculate_GreenR2, except for sum of elements
  !Note: to backup old version zconcat calls (and Glsub deallocations) must be
  !      uncommented and all this part removed
  !If only one block is present, concatenation is not needed and it's implemented in a
  !more trivial way
  subroutine blk2csr(G,struct,P,Gcsr)

    type(z_DNS), dimension(:,:) :: G
    type(Tstruct_info), intent(in) :: struct
    type(z_CSR) :: Gcsr
    type(z_CSR) :: P, G_sp

    integer :: nbl, oldx, row, col, iy, ix, x, y, ii, jj, nrows

    nbl = struct%num_PLs
    nrows = struct%mat_PL_end(nbl)

    !create Gcsr with same pattern of P
    call create(Gcsr,P%nrow,P%ncol,P%nnz)
    Gcsr%rowpnt = P%rowpnt
    Gcsr%colind = P%colind
    Gcsr%nzval = (0.0_dp, 0.0_dp)

    associate(indblk=>struct%mat_PL_start)
    !Cycle upon all rows
    x = 1
    do ii = 1, nrows
      !Search block x containing row ii
      oldx = x
      if (oldx.EQ.nbl) THEN
        x = oldx
      ELSE
        do ix = oldx, oldx+1
          if ( (ii.GE.indblk(ix)).AND.(ii.LT.indblk(ix+1)) ) x = ix
        end do
      endif

      !Offset: row is the index for separate blocks
      row = ii - indblk(x) + 1

      !Cycle upon columns of Gcsr (which has been ALREADY MASKED by S)
      do jj = Gcsr%rowpnt(ii), Gcsr%rowpnt(ii+1) -1
        if (Gcsr%colind(jj).gt.nrows) CYCLE
        !Choose which block column we're dealing with
        y = 0
        if (x.eq.1) then
          if ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then
            y = 1
          ELSEif ( (Gcsr%colind(jj).GE.indblk(x + 1)).AND.(Gcsr%colind(jj).LT.indblk(x + 2)) ) then
            y = 2
          endif
        elseif (x.eq.nbl) then
          if ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then
            y = nbl
          ELSEif ( (Gcsr%colind(jj).GE.indblk(x - 1)).AND.(Gcsr%colind(jj).LT.indblk(x)) ) then
            y = nbl - 1
          endif
        else
          do iy = x-1, x+1
            if ( (Gcsr%colind(jj).GE.indblk(iy)).AND.(Gcsr%colind(jj).LT.indblk(iy + 1)) ) y = iy
          end do
        endif

        if (y.EQ.0) THEN
          write(*,*)
          write(*,*) 'ERROR in blk2csr: probably wrong PL size', x
          write(*,*) 'row',ii,nrows,Gcsr%colind(jj)
          write(*,*) 'block indeces:',indblk(1:nbl)
          STOP
        endif

        col = Gcsr%colind(jj) - indblk(y) + 1

        if (allocated(G(x,y)%val)) THEN
          Gcsr%nzval(jj) = Gcsr%nzval(jj) + G(x,y)%val(row,col)
        endif

      end do

    end do
    end associate

  end subroutine blk2csr


  !****************************************************************************
  !
  !  Calculate Green Retarded in the
  !  contacts regions, where overlap with device orbitals is non-zero
  !  (only the upper part, needed in both K=0 and K points calculations,
  !   or both upper and lower parts)
  !  writing on memory
  !
  !****************************************************************************

  subroutine calculate_Gr_outer(Tlc,Tcl,gsurfR,struct,lower,Aout)

    !****************************************************************************
    !Input:
    !Tlc: sparse arraymatrices array containing contact-device interacting blocks
    !     ESH-H
    !Tlc: sparse arraymatrices array containing device-contact interacting blocks
    !     ESH-H
    !gsurfR: sparse matrices array containing contacts surface green
    !lower: if .true., also lower parts are calculated and concatenated
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont),
    !cindblk(ncont), Gr in the diagonal block related to the interaction layer
    !
    !Output:
    !Aout: sparse matrix containing density matrix in the region
    !      corresponding to non-zero overlap
    !
    !****************************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:) :: Tlc,Tcl,gsurfR
    logical :: lower
    type(Tstruct_info), intent(in) :: struct
    type(z_CSR) :: Aout


    !Work
    type(z_DNS) :: work1, Grcl, Grlc
    type(z_CSR) :: GrCSR, TCSR
    integer :: i,cb,nrow_tot,i1,j1
    integer :: ncont, nbl

    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim

    if (.not.allocated(Aout%nzval)) THEN
      call create(Aout,nrow_tot,nrow_tot,0)
      Aout%rowpnt(:)=1
    endif

    do i=1,ncont

      !Numero di blocco del contatto
      cb=struct%cblk(i)
      call prealloc_mult(Gr(cb,cb),Tlc(i),(-1.0_dp, 0.0_dp),work1)
      call prealloc_mult(work1,gsurfR(i),Grlc)

      call destroy(work1)

      j1=nzdrop(Grlc,EPS)
      call create(GrCSR,Grlc%nrow,Grlc%ncol,j1)
      call dns2csr(Grlc,GrCSR)
      call destroy(Grlc)
      j1=nzdrop(Tlc(i),EPS)
      call create(TCSR,Tlc(i)%nrow,Tlc(i)%ncol,j1)
      call dns2csr(Tlc(i),TCSR)
      call zmask_realloc(GrCSR,TCSR)
      call destroy(TCSR)

      !Concatenazione di Asub nella posizione corrispondente
      i1=struct%mat_PL_start(cb)
      j1=struct%mat_B_start(i)

      call concat(Aout,GrCSR,i1,j1)

      call destroy(GrCSR)

      if (lower) THEN

        call prealloc_mult(gsurfR(i),Tcl(i),(-1.0_dp, 0.0_dp), work1)
        call prealloc_mult(work1, Gr(cb,cb), Grcl)

        call destroy(work1)

        j1=nzdrop(Grcl,EPS)
        call create(GrCSR,Grcl%nrow,Grcl%ncol,j1)
        call dns2csr(Grcl,GrCSR)
        call destroy(Grcl)
        j1=nzdrop(Tcl(i),EPS)
        call create(TCSR,Tcl(i)%nrow,Tcl(i)%ncol,j1)
        call dns2csr(Tcl(i),TCSR)
        call zmask_realloc(GrCSR,TCSR)
        call destroy(TCSR)

        i1 = struct%mat_B_start(i)-struct%central_dim+struct%mat_PL_start(nbl+1)-1
        j1 = struct%mat_PL_start(cb)

        call concat(Aout,GrCSR,i1,j1)

        call destroy(GrCSR)

      endif

    end do

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'Outer_GreenR done'
      WRITE(*,*) '********************'
    endif

  end subroutine calculate_Gr_outer



  !****************************************************************************
  !
  !  Calculate Gn contributions for all contacts except collector, in the
  !  outer region where contacts-device overlap is non-zero - writing on
  !  memory
  !
  !****************************************************************************

  subroutine calculate_Gn_outer(Tlc,gsurfR,SelfEneR,struct,frm,ref,lower,Glout)

    !****************************************************************************
    !Input:
    !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
    !gsurfR: sparse matrices array containing contacts surface green
    !SelfEneR: sparse matrices array containing contacts Self Energy
    !frm: array containing Fermi distribution values for all contacts
    !ref: reference contact
    !
    !global variables needed: nbl, indblk(nbl+1), cindblk(ncont), ncont,
    !cblk(ncont), Gr(:,:), diagonal, subdiagonal, overdiagonal and
    !in colums Gr(:,cb) where cb are layers interacting with all contacts
    !but collector
    !
    !Output:
    !Glout: sparse matrix containing G_n  contributions in the region
    !       corresponding to non-zero overlap
    !
    !****************************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:) :: Tlc, gsurfR, SelfEneR
    real(dp), dimension(:) :: frm
    type(Tstruct_info), intent(in) :: struct
    integer :: ref
    logical :: lower
    type(z_CSR) :: Glout

    !Work
    type(z_DNS) :: Gam, gsurfA, Ga, work1, work2, work3, Glsub
    type(z_CSR) :: GlCSR, TCSR
    integer :: j,k,cb,cbj,i1,j1,nrow_tot
    integer :: ncont, nbl
    complex(dp) :: frmdiff


    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim

    !Allocazione della Glout
    !Righe totali del conduttore effettivo nrow_tot
    !nrow_tot=indblk(nbl+1)-1
    !do i=1,ncont
    !   nrow_tot=nrow_tot+ncdim(i) !gsurfR(i)%nrow
    !end do
    if (.not.allocated(Glout%nzval)) THEN
      call create(Glout,nrow_tot,nrow_tot,0)
      Glout%rowpnt(:)=1
    endif
    !***
    !Iterazione su tutti i contatti "k"
    !***
    do k=1,ncont

      !Esegue le operazioni relative al contatto solo se e` valida la condizione
      !sulle distribuzioni di Fermi e se non si tratta del contatto iniettante (ref)
      if ((ABS(frm(k)-frm(ref)).GT.EPS).AND.(k.NE.ref)) THEN

        !Calcolo della Gamma corrispondente
        cb=struct%cblk(k)
        !nota:cb indica l'indice di blocco corrispondente al contatto j-esimo

        call zspectral(SelfEneR(k),SelfEneR(k),0,Gam)

        !***
        !Calcolo del contributo sulla proria regione
        !***
        frmdiff = cmplx(frm(ref)-frm(k),0.0_dp,dp)

        call zspectral(gsurfR(k),gsurfR(k),0,work1)

        !print*, 'work2=Tlc*work1=Tlc*j(gsurfR-gsurfA)'
        call prealloc_mult(Tlc(k),work1,work2)
        call destroy(work1)

        !print *, 'work1=-Gr(cb,cb)*work2=-Gr(cb,cb)*Tlc*j(gsurfR-gsurfA)'
        call prealloc_mult(Gr(cb,cb),work2,frmdiff,work1)
        call destroy(work2)

        !print*,'work2=Tlc*gsurfA'
        call zdagger(gsurfR(k),gsurfA)
        call prealloc_mult(Tlc(k),gsurfA,work2)
        call destroy(gsurfA)

        !print*,'work3=Ga*work2=Ga*Tlc*gsurfA'
        call zdagger(Gr(cb,cb),Ga)
        call prealloc_mult(Ga,work2,work3)

        call destroy(Ga)
        call destroy(work2)

        !print*,'work2=Gam*work3=Gam*Ga*Tlc*gsurfA'
        call prealloc_mult(Gam,work3,work2)
        call destroy(work3)

        !print *,'work3=-Gr*work2=-Gr*Gam*Ga*Tlc*gsurfA'
        call prealloc_mult(Gr(cb,cb),work2,frmdiff,work3)
        call destroy(work2)

        !Contributo totale sulla propria regione
        call prealloc_sum(work3,work1,Glsub)
        call destroy(work1)
        call destroy(work3)

        call mask(Glsub,Tlc(k))
        i1=nzdrop(Glsub,EPS)

        if (i1.gt.0) THEN
          call create(GlCSR,Glsub%nrow,Glsub%ncol,i1)
          call dns2csr(Glsub,GlCSR)

          !Concatenazione di Glsub nella matrice globale Glout
          i1=struct%mat_PL_start(cb)
          j1=struct%mat_B_start(k)-struct%central_dim+struct%mat_PL_start(nbl+1)-1
          call concat(Glout,GlCSR,i1,j1)

          ! compute lower outer part using (iG<)+ = iG<
          if (lower) THEN
            call zdagger(GlCSR,TCSR)
            call concat(Glout,TCSR,j1,i1)
            call destroy(TCSR)
          endif

          call destroy(GlCSR)
        end if

        call destroy(Glsub)
        !***
        !Ciclo per il calcolo dei contributi sulle regioni degli altri contatti
        !***

        do j=1,ncont

          cbj=struct%cblk(j)
          !Esegue le operazioni del ciclo solo se il j.ne.k o se
          !il blocco colonna di Gr e` non nullo (altrimenti il contributo e` nullo)

          if ((j.NE.k).AND.(Gr(cbj,cb)%nrow.NE.0 .AND. (Gr(cbj,cb)%ncol.NE.0))) THEN

            !print*,'work1=Tlc*gsurfA'
            call zdagger(gsurfR(j),gsurfA)
            call prealloc_mult(Tlc(j),gsurfA,work1)
            call destroy(gsurfA)

            !print*,'work2=Ga*work1=Ga*Tlc*gsurfA'
            call zdagger(Gr(cbj,cb),Ga)
            call prealloc_mult(Ga,work1,work2)

            call destroy(Ga)
            call destroy(work1)

            !print*,'work1=Gam*work2=Gam*Ga*Tlc*gsurfA'
            call prealloc_mult(Gam,work2,work1)
            call destroy(work2)

            !print*,'Glsub=-Gr*work1=-Gr*Gam*Ga*Tlc*gsurfA'
            call prealloc_mult(Gr(cbj,cb),work1,frmdiff,Glsub)
            call destroy(work1)

            call mask(Glsub,Tlc(j))
            i1=nzdrop(Glsub,EPS)

            if (i1.gt.0) THEN
              call create(GlCSR,Glsub%nrow,Glsub%ncol,i1)
              call dns2csr(Glsub,GlCSR)

              !Concatenazione di Glsub nella posizione corrispondente al contatto "j"
              i1=struct%mat_PL_start(cbj)
              j1=struct%mat_B_start(j)-struct%central_dim+struct%mat_PL_start(nbl+1)-1

              call concat(Glout,GlCSR,i1,j1)

              ! compute lower outer part using (iG<)+ = iG<
              if (lower) THEN
                call zdagger(GlCSR,TCSR)
                call concat(Glout,TCSR,j1,i1)
                call destroy(TCSR)
              endif

              call destroy(GlCSR)
            endif

            call destroy(Glsub)

          endif
        end do

        call destroy(Gam)

      endif

    end do

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'Outer_Gl done'
      WRITE(*,*) '********************'
    endif

  end subroutine calculate_Gn_outer


  !---------------------------------------------------

  subroutine calculate_transmissions(H,S,Ec,SelfEneR,ni,nf,str,tun_proj,tun_mat)
    Type(z_CSR) :: H
    Type(z_CSR) :: S
    Complex(dp) :: Ec
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    Integer :: ni(:)
    Integer :: nf(:)
    Type(z_CSR) :: ESH_tot
    Type(TStruct_Info) :: str
    type(intarray), intent(in) :: tun_proj
    Real(dp), Dimension(:) :: tun_mat

    ! Local variables
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Real(dp) :: tun
    Integer :: nbl,ncont,ibl
    Integer :: i, ierr, icpl, nit, nft, nt, nt1

    nbl = str%num_PLs
    ncont = str%num_conts

    if (ncont == 1) then
      tun_mat = 0.0_dp
      return
    end if

    !Calculation of ES-H and brak into blocks
    call prealloc_sum(H,S,(-1.0_dp, 0.0_dp),Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
      ibl = str%cblk(i)
      ESH(ibl,ibl)%val = ESH(ibl,ibl)%val-SelfEneR(i)%val
    end do

    nit=ni(1)
    nft=nf(1)
    ! find the contact with smaller block index
    if (str%cblk(nit).lt.str%cblk(nft)) then
      nt = str%cblk(nit)
    else
      nt = str%cblk(nft)
    endif

    ! Fall here when there are 2 contacts for fast transmission
    if (ncont == 2 .and. size(ni) == 1 .and. nt == 1) then
      call allocate_gsm(gsmr,nbl)
      call calculate_gsmr_blocks(ESH,nbl,2,.false.)
      call allocate_blk_dns(Gr,nbl)
      call calculate_Gr_tridiag_blocks(ESH,1)
      call calculate_single_transmission_2_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun_proj,tun)
      tun_mat(1) = tun

    else
      ! MULTITERMINAL case
      call allocate_gsm(gsmr,nbl)
      call calculate_gsmr_blocks(ESH,nbl,2)

      do icpl = 1, size(ni)

        !Computation of transmission(s) between contacts ni(:) -> nf(:)
        nit=ni(icpl)
        nft=nf(icpl)

        ! find the largest contact block between the two terminals
        if (str%cblk(nit).gt.str%cblk(nft)) then
          nt1 = str%cblk(nit)
        else
          nt1 = str%cblk(nft)
        endif

        if (icpl == 1) then
          ! Iterative calculation of Gr down to nt
          nt = nt1
          call allocate_blk_dns(Gr,nbl)
          call calculate_Gr_tridiag_blocks(ESH,1)
          if (nt.gt.1) then
            call calculate_Gr_tridiag_blocks(ESH,2,nt)
          end if
        else
          ! When more contacts are present sometimes we can re-use previous GF
          ! if nt1 > nt extend the Gr calculation
          if (nt1 .gt. nt) then
            call calculate_Gr_tridiag_blocks(ESH,nt+1,nt1)
            nt = nt1
          endif
        end if

        call calculate_single_transmission_N_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun_proj,tun)

        tun_mat(icpl) = tun

      end do
    end if

    !Distruzione delle Green
    do i=2, nt
      call destroy(Gr(i,i))
      call destroy(Gr(i-1,i))
      call destroy(Gr(i,i-1))
    end do
    call destroy(Gr(1,1))

    call deallocate_blk_dns(Gr)

    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

  end subroutine calculate_transmissions

  !************************************************************************
  !
  ! Subroutine for transmission calculation
  !
  !************************************************************************
  ! NOTE:
  !
  !  This subroutine was hacked quickly to obain effiecient tunneling calcs
  !  Useful only when there are 2 contacts
  !                ===================
  !************************************************************************

  subroutine calculate_single_transmission_2_contacts(ni,nf,ESH,SelfEneR,cblk,tun_proj,TUN)
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    real(dp), intent(out) :: TUN

    !Work variables
    Integer :: ct1, bl1
    logical, dimension(:), allocatable :: tun_mask
    Type(z_DNS) :: work1, work2, GAM1_dns, GA, TRS, AA
    Complex(dp), parameter ::    j = (0.0_dp,1.0_dp)  ! CMPX unity

    if (size(cblk).gt.2) then
      write(*,*) "ERROR: calculate_single_transmission_2_contacts is valid only for 2 contacts"
      TUN = 0.0_dp
      return
    endif

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
      ct1=ni;
    else
      ct1=nf;
    endif

    bl1=cblk(ct1);

    call zdagger(Gr(bl1,bl1),GA)

    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call prealloc_mult(GAM1_dns,Gr(bl1,bl1),work1)

    call prealloc_mult(work1,GAM1_dns,work2)

    call destroy(work1)

    call prealloc_mult(work2,GA,work1)

    call destroy(work2)

    call create(AA,GA%nrow,GA%ncol)

    AA%val = j * (Gr(bl1,bl1)%val-GA%val)

    call destroy(GA)

    call prealloc_mult(GAM1_dns,AA,work2)

    call destroy(GAM1_dns,AA)

    call create(TRS,work1%nrow,work1%ncol)

    TRS%val = work2%val - work1%val

    call get_tun_mask(ESH, bl1, tun_proj, tun_mask)

    TUN = abs( real(trace(TRS, tun_mask)) )

    call log_deallocate(tun_mask)

    call destroy(TRS,work1,work2)

  end subroutine calculate_single_transmission_2_contacts

  !************************************************************************
  !
  ! Subroutine for transmission calculation (GENERIC FOR N CONTACTS)
  !
  !************************************************************************
  subroutine calculate_single_transmission_N_contacts(ni,nf,ESH,SelfEneR,cblk,tun_proj,TUN)
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    real(dp), intent(out) :: TUN

    !Work variables
    Integer :: ct1, ct2, bl1, bl2, i, nbl
    logical, dimension(:), allocatable :: tun_mask
    Type(z_DNS) :: work1, work2, GAM1_dns, GAM2_dns, GA, TRS
    Real(kind=dp) :: max

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
      ct1=ni;ct2=nf;
    else
      ct1=nf;ct2=ni;
    endif

    bl1=cblk(ct1); bl2=cblk(ct2);
    nbl = size(cblk)
    ! in this way nt1 < nt2 by construction

    if ( nbl.gt.1 .and. (bl2-bl1).gt.1) then

      ! Compute column-blocks of Gr(i,bl1) up to i=bl2
      ! Gr(i,bl1) = -gr(i) T(i,i-1) Gr(i-1,bl1)
      do i = bl1+1, bl2
        !Checks whether previous block is non null.
        !If so next block is also null => TUN = 0
        max=maxval(abs(Gr(i-1,bl1)%val))

        if (max.lt.EPS) then
          TUN = EPS*EPS !for log plots
          !Destroy also the block adjecent to diagonal since
          !this is not deallocated anymore in calling subroutine
          if (i.gt.(bl1+1)) call destroy(Gr(i-1,bl1))
          return
        endif

        !Checks whether block has been created, if not do it
        if (.not.allocated(Gr(i,bl1)%val)) then

          call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.0_dp, 0.0_dp),work1)

          call prealloc_mult(work1,Gr(i-1,bl1),Gr(i,bl1))

          call destroy(work1)

        endif

        ! avoid destroying blocks closer to diagonal
        if (i.gt.(bl1+2)) call destroy(Gr(i-1,bl1))

      end do

    endif

    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)
    call zspectral(SelfEneR(ct2),SelfEneR(ct2),0,GAM2_dns)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call prealloc_mult(GAM2_dns,Gr(bl2,bl1),work1)

    call destroy(GAM2_dns)

    call prealloc_mult(work1,GAM1_dns,work2)

    call destroy(work1)

    call destroy(GAM1_dns)

    call zdagger(Gr(bl2,bl1),GA)

    if (bl2.gt.bl1+1) call destroy( Gr(bl2,bl1) )

    call prealloc_mult(work2,GA,TRS)

    call destroy(work2)

    call destroy(GA)

    call get_tun_mask(ESH, bl2, tun_proj, tun_mask)

    TUN = abs( real(trace(TRS, tun_mask)) )

    call log_deallocate(tun_mask)

    call destroy(TRS)

  end subroutine calculate_single_transmission_N_contacts

  ! Based on projection indices build a logical mask just on contact block
  subroutine get_tun_mask(ESH,nbl,tun_proj,tun_mask)
    Type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: nbl
    type(intarray), intent(in) :: tun_proj
    logical, intent(out), allocatable :: tun_mask(:)

    integer :: ii, istart, iend, ind

    call log_allocate(tun_mask, ESH(nbl,nbl)%nrow)

    if (allocated(tun_proj%indexes)) then
      tun_mask = .false.

      ! set the start/end indices of nbl
      ! NB: istart has offset -1 to avoid +/-1 operations
      istart = 0
      do ii = 1, nbl-1
        istart = istart + ESH(ii,ii)%nrow
      end do
      iend = istart + ESH(nbl,nbl)%nrow + 1

      ! select the indices in tun_proj
      do ii = 1, size(tun_proj%indexes)
         ind = tun_proj%indexes(ii)
         if (ind > istart .and. ind < iend) then
            tun_mask(ind - istart) = .true.
         end if
      end do
    else
      tun_mask = .true.
    end if

  end subroutine get_tun_mask

  !---------------------------------------------------!
  !Subroutine for transmission and dos calculation    !
  !---------------------------------------------------!

  subroutine calculate_transmissions_and_dos(H,S,Ec,SelfEneR,Gs,ni,nf,str,tun_proj,TUN_MAT,dos_proj,LEDOS)
    Type(z_CSR), intent(in) :: H
    Type(z_CSR), intent(in) :: S
    Complex(dp), intent(in) :: Ec
    Type(z_DNS), Dimension(MAXNCONT), intent(in) :: SelfEneR, Gs
    Integer, intent(in) :: ni(:)
    Integer, intent(in) :: nf(:)
    Type(TStruct_Info), intent(in) :: str
    type(intarray), intent(in) :: tun_proj
    type(intarray), dimension(:), intent(in) :: dos_proj
    Real(dp), Dimension(:), intent(inout) :: TUN_MAT
    Real(dp), Dimension(:), intent(inout) :: LEDOS

    ! Local variables
    Type(z_CSR) :: ESH_tot, GrCSR
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Type(r_CSR) :: Grm                          ! Green Retarded nella molecola
    real(dp), dimension(:), allocatable :: diag
    Real(dp) :: tun
    Complex(dp) :: zc
    Integer :: nbl,ncont, ierr
    Integer :: nit, nft, icpl
    Integer :: iLDOS, i2, i
    Character(1) :: Im


    nbl = str%num_PLs
    ncont = str%num_conts
    Im = 'I'

    !Calculation of ES-H and brak into blocks
    call prealloc_sum(H,S,(-1.0_dp, 0.0_dp),Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
      ESH(str%cblk(i),str%cblk(i))%val = ESH(str%cblk(i),str%cblk(i))%val-SelfEneR(i)%val
    end do

    call allocate_gsm(gsmr,nbl)
    call calculate_gsmr_blocks(ESH,nbl,2)

    call allocate_blk_dns(Gr,nbl)
    call calculate_Gr_tridiag_blocks(ESH,1)
    call calculate_Gr_tridiag_blocks(ESH,2,nbl)

    !Computation of transmission(s) between contacts ni(:) -> nf(:)
    do icpl=1,size(ni)

      nit=ni(icpl)
      nft=nf(icpl)

      select case(ncont)
      case(1)
        tun = 0.0_dp
      case(2)
        call calculate_single_transmission_2_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun_proj,tun)
      case default
        call calculate_single_transmission_N_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun_proj,tun)
      end select

      TUN_MAT(icpl) = tun

    end do

    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    !Distruzione dei blocchi fuori-diagonale
    do i=2,nbl
      call destroy(Gr(i-1,i))
      call destroy(Gr(i,i-1))
    end do
    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

    call create(Grm,str%mat_PL_start(nbl+1)-1,str%mat_PL_start(nbl+1)-1,0)

    Grm%rowpnt(:)=1

    do i=1,nbl
      call create(GrCSR,Gr(i,i)%nrow,Gr(i,i)%ncol,Gr(i,i)%nrow*Gr(i,i)%ncol)
      call dns2csr(Gr(i,i),GrCSR)
      !Concatena direttamente la parte immaginaria per il calcolo della doS
      zc=(-1.0_dp,0.0_dp)/pi

      call concat(Grm,zc,GrCSR,Im,str%mat_PL_start(i),str%mat_PL_start(i))
      call destroy(Gr(i,i))
      call destroy(GrCSR)
    end do

    call deallocate_blk_dns(Gr)

    !Compute LDOS on the specified intervals
    if (size(dos_proj).gt.0) then
      call log_allocate(diag, Grm%nrow)
      call getdiag(Grm,diag)
      do iLDOS=1,size(dos_proj)
        do i = 1, size(dos_proj(iLDOS)%indexes)
          i2 = dos_proj(iLDOS)%indexes(i)
          if (i2 .le. str%central_dim) then
            LEDOS(iLDOS) = LEDOS(iLDOS) + diag(i2)
          end if
        end do
      end do
      call log_deallocate(diag)
    endif

    call destroy(Grm)

  end subroutine calculate_transmissions_and_dos


  !---------------------------------------------------


  subroutine allocate_gsm(gsm,nbl)
    type(z_DNS), dimension(:), allocatable :: gsm
    integer :: nbl, ierr

    allocate(gsm(nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate gsm'

  end subroutine allocate_gsm

  !---------------------------------------------------

  subroutine allocate_blk_dns(blkM,nbl)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: nbl, ierr

    allocate(blkM(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'

  end subroutine allocate_blk_dns

  !---------------------------------------------------

  subroutine deallocate_gsm(gsm)
    type(z_DNS), dimension(:), allocatable :: gsm
    integer :: ierr

    deallocate(gsm,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate gsmr'

  end subroutine deallocate_gsm

  !---------------------------------------------------

  subroutine deallocate_blk_dns(blkM)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: ierr

    deallocate(blkM,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate block-Matrix'

  end subroutine deallocate_blk_dns

end module iterative
