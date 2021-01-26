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

! Inelastic electron-phonon interactions
! This is the base class for inelastic interactions

module elphinel

  use ln_precision, only : dp
  use interactions, only : interaction
  use ln_allocation, only : log_allocate, log_deallocate
  use ln_structure, only : TStruct_info
  use ln_structure, only : TBasisCenters, create_TBasis, destroy_TBasis
  use mat_def, only : z_csr, z_dns, create, destroy
  use ln_cache
  use distributions, only : bose
  use iso_c_binding, only : c_int, c_double

  implicit none
  private

  public :: ElPhonInel, ElPhonInel_create

  type, extends(interaction) :: ElPhonInel
    private
    !> holds atomic structure
    type(TBasisCenters) :: basis

    !> Bose-Einstein phonon occupation
    real(dp) :: Nq

    !> sigma_r and sigma_n
    type(TMatrixCache), allocatable :: sigma_r
    type(TMatrixCache), allocatable :: sigma_n
    !> Gr and Gn are pointer alias stored in negf.
    type(TMatrixCache), pointer :: G_r => null()
    type(TMatrixCache), pointer :: G_n => null()

  contains

    procedure :: add_sigma_r
    procedure :: get_sigma_n
    procedure :: set_Gr
    procedure :: set_Gn
    procedure :: compute_sigma_r
    procedure :: compute_sigma_n

  end type ElPhonInel

  !> Interface of C function that perform all the MPI communications in order to compute
  !  sigma_r and sigma_n
  interface
   integer (c_int) function self_energy(Np, NK, NE, NKloc, NEloc, iEshift, GG, sigma, &
            & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, cart_comm, fac_min, fac_plus, func) &
            & bind(C, name='self_energy')
     import :: dp
     !> matrix size Np x Np
     integer(c_int), value :: Np
     !> Total Number of K-points
     integer(c_int), value :: NK
     !> Total Number of E-points
     integer(c_int), value :: NE
     !> Total Number of local K-points
     integer(c_int), value :: NKloc
     !> Total Number of local E-points
     integer(c_int), value :: NEloc
     !> grid shift of hwq  iEshif
     integer(c_int), value :: iEshift
     !> Green Function
     complex(c_double_complex) :: GG
     !> self energy
     complex(c_double_complex) :: Self energy
     !> 6 Buffers of size Np x Np
     complex(c_double_complex) :: sbuff1
     complex(c_double_complex) :: sbuff2
     complex(c_double_complex) :: rbuff1
     complex(c_double_complex) :: rbuff2
     complex(c_double_complex) :: sbuffH
     complex(c_double_complex) :: rbuffH
     !> Cartesian MPI_communicator (negf%cartgrid%id)
     integer(c_int) :: cart_comm
     !> prefactor of GG(E-wq) (e.g. nq+1, nq, i/2)
     complex(c_double_complex) :: fac_min
     !> prefactor of GG(E+wq)
     complex(c_double_complex) :: fac_plus
     !> call-back function K(k,q) for integration
     type(C_FUNPTR), value :: callback
   end subroutine self_energy

contains

  !>
  ! Factory for el-ph inelastic model based on polar-optical modes
  ! @param struct : contact/device partitioning infos
  ! @param coupling: coupling (energy units)
  ! @param wq: phonon frequence
  ! @param niter: fixed number of scba iterations
  subroutine ElPhonInel_create(this, struct, basis, coupling, wq, Temp, niter)

    type(ElPhonInel), intent(inout) :: this
    type(TStruct_Info), intent(in) :: struct
    type(TBasisCenters), intent(in) :: basis
    real(dp), dimension(:), intent(in) :: coupling
    real(dp), intent(in) :: wq
    real(dp), intent(in) :: Temp
    integer, intent(in) :: niter

    this%descriptor = &
        & "Electron-Phonon inelastic model for optical phonons"

    this%scba_niter = niter
    this%struct = struct
    this%basis = basis
    this%wq = wq
    this%Nq = bose(wq, Temp)

    ! Initialize the cache space
    this%sigma_r = TMatrixCacheMem()
    this%sigma_n = TMatrixCacheMem()

  end subroutine ElPhonInel_create


  !> This interface should append
  !  the retarded self energy to ESH
  subroutine add_sigma_r(this, esh, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: n, npl, ii, ierr, jj
    type(z_dns) :: tmp_blk
    type(TMatLabel) :: label

    if (this%scba_iter .eq. 0) return

    ! This could be done more performant, but this model won't see much use
    ! so I am leaving this way
    label%kpoint = k_index
    label%energy_point = en_index
    label%spin = 0

    do jj = 1, npl
      label%row_block = jj
      label%col_block = jj
      call sigma_r%retrieve(tmp_blk, label)
      ESH(jj, jj)%val = ESH(jj, jj)%val - tmp_blk%val
      call destroy(tmp_blk)
      if (jj .lt. npl) then
        label%row_block = jj
        label%col_block = jj + 1
        call sigma_r%retrieve(tmp_blk, label)
        ESH(jj, jj + 1)%val = ESH(jj, jj + 1)%val - tmp_blk
        call destroy(tmp_blk)
        ESH(jj + 1, jj)%val = ESH(jj + 1, jj)%val - tmp_blk
        call destroy(tmp_blk)
      end if
    end do

  end subroutine add_sigma_r


  !> Returns the lesser (n) Self Energy in block format
  !
  subroutine get_sigma_n(this, blk_sigma_n, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    if (this%scba_iter .eq. 0) return

  end subroutine get_sigma_n


  !> Give the Gr at given energy point to the interaction
  subroutine set_Gr(this, Gr, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    ! here we need to bind the pointer in negf


  end subroutine set_Gr


  !> Give the Gn at given energy point to the interaction
  subroutine set_Gn(this, Gn, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    ! here we need to bind the pointer in negf

  end subroutine set_Gn


  !> Give the Gn at given energy point to the interaction
  subroutine compute_Sigma_r(this, en_index, k_index, spin)
    class(ElPhonInel) :: this
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: cart_comm
    integer :: ibl, nbl, Np, NK, NE, NKloc, NEloc, iEshift
    real(c_double) :: factor_min, factor_plus
    complex(c_double_complex), allocatable :: sbuff1(:,:), rbuff1(:,:)
    complex(c_double_complex), allocatable :: sbuff2(:,:), rbuff2(:,:)
    complex(c_double_complex), allocatable :: sbuffH(:,:), rbuffH(:,:)

    nbl = this%struct%num_PLs

    do ibl = 1, nbl

      ! block dimension
      Np = this%struct%mat_PL_start(ibl) - this%struct%mat_PL_end(ibl) + 1

      ! create buffers
      call log_allocate(sbuff1, Np, Np)
      call log_allocate(rbuff1, Np, Np)
      call log_allocate(sbuff2, Np, Np)
      call log_allocate(rbuff2, Np, Np)
      call log_allocate(sbuffH, Np, Np)
      call log_allocate(rbuffH, Np, Np)

      NK =
      NKloc =
      NE =
      NEloc =

      iEshift = nint(this%wq/dE)

      factor_minus = this%Nq + 1
      factor_plus = this%Nq

      mat_Gr = G_r%retrieve(iE, ik, ispin, ibl, ibl)

      ! Call to SEBASTIAN c-function to get GG
      self_energy(Np, NK, NE, NKloc, NEloc, iEshift, mat_Gr, sigma_r, &
            & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, &
            & cart_comm, factor_min, factor_plus, kappa)

    end do

  end subroutine compute_Sigma_r


  !> Give the Gn at given energy point to the interaction
  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(ElPhonInel) :: this
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    ! here we need to bind the pointer in negf


  end subroutine compute_Sigma_n


  !> Integral over qz of the electron-phonon polar optical couping.
  real(c_double) function kappa() bind(c, name='kappa')



  end function kappa


end module elphinel
