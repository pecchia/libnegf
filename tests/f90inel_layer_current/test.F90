!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Green's functions.        !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Function calculation                             !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute it and/or modify          !
!! it under the terms of the GNU Lesser General Public License as published !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!

!> Inelastic SCBA with polar-optical phonons on the 100-site chain of f90_transmission.
!>
!> Layer currents and Meir-Wingreen currents are computed on the same energy grid twice:
!> on all tasks with one energy point per task, and on one task holding all the points.
!> The two must agree.
program test

  use libnegf
  use lib_param
  use integrations
  use mpi_f08

  implicit none

  integer, parameter :: wp = kind(1.d0)
  !> Relative tolerance between the runs on all tasks and on one task
  real(wp), parameter :: tolerance = 1.0e-10_wp
  !> Relative tolerances of current conservation and of the comparison with the reference
  real(wp), parameter :: conservationTol = 1.0e-5_wp, referenceTol = 1.0e-6_wp

  !> Reference currents of the run with 4 tasks, hence 4 energy points
  real(wp), parameter :: refLayerCurrents(5) = [-3.2971874380624597_wp, -3.2971874845242626_wp, &
      & -3.2971875108932482_wp, -3.2971875165512410_wp, -3.2971875026442721_wp]
  real(wp), parameter :: refMeirWingreen(2) = [-3.2971875426524249_wp, 3.2971873447110109_wp]

  type(MPI_Comm) :: singleComm
  integer :: ierr, rank, nRanks, color

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nRanks, ierr)

  if (nRanks < 2) then
    call fail('This test needs at least 2 MPI tasks')
  end if

  color = merge(0, MPI_UNDEFINED, rank == 0)
  call MPI_Comm_split(MPI_COMM_WORLD, color, rank, singleComm, ierr)

  call check_distributed_against_serial(.true.)
  call check_distributed_against_serial(.false.)

  if (rank == 0) write(*,*) 'Done'
  call MPI_Finalize(ierr)

contains

  !> Run on all tasks (one point per task) and on one task (all points), and compare
  subroutine check_distributed_against_serial(tLayerCurrents)
    logical, intent(in) :: tLayerCurrents

    real(wp), allocatable :: distributed(:,:), serial(:,:), total(:)
    real(wp) :: scale
    character(len=32) :: name

    name = merge('layer currents', 'Meir-Wingreen ', tLayerCurrents)

    call compute_currents(MPI_COMM_WORLD, nRanks, tLayerCurrents, distributed)
    if (rank == 0) then
      call compute_currents(singleComm, nRanks, tLayerCurrents, serial)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (rank /= 0) return

    write(*,*) trim(name), ': current per interface (summed over energies)'
    total = sum(distributed, dim=1)
    scale = maxval(abs(total))
    write(*,'(10ES25.16)') total
    write(*,*) trim(name), ': max difference distributed/serial ', maxval(abs(distributed - serial))

    if (scale < 1.0e-8_wp) then
      call fail('Vanishing current in '//trim(name))
    end if
    if (maxval(abs(distributed - serial)) > tolerance * maxval(abs(distributed))) then
      call fail('Distributed and serial runs give different '//trim(name))
    end if

    ! The current is the same at every interface, and opposite at the two contacts
    if (tLayerCurrents) then
      if (maxval(total) - minval(total) > conservationTol * scale) then
        call fail('Layer currents are not conserved')
      end if
    else
      if (abs(sum(total)) > conservationTol * scale) then
        call fail('Meir-Wingreen currents are not conserved')
      end if
    end if

    if (nRanks == 4) then
      if (tLayerCurrents) then
        if (maxval(abs(total - refLayerCurrents)) > referenceTol * scale) then
          call fail('Layer currents differ from the reference')
        end if
      else
        if (maxval(abs(total - refMeirWingreen)) > referenceTol * scale) then
          call fail('Meir-Wingreen currents differ from the reference')
        end if
      end if
    end if

  end subroutine check_distributed_against_serial


  !> Compute the currents with the inelastic SCBA on an energy grid of nPoints points
  subroutine compute_currents(comm, nPoints, tLayerCurrents, currents)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: nPoints
    logical, intent(in) :: tLayerCurrents
    !> Current for each energy and interface (layer currents) or contact (Meir-Wingreen)
    real(wp), allocatable, intent(out) :: currents(:,:)

    integer, parameter :: nAtoms = 100, nCentral = 60
    real(wp), parameter :: siteSpacing = 0.5_wp, Estep = 0.05_wp, kT = 1.0e-3_wp
    real(wp), parameter :: Emin = -0.1_wp

    ! Phonon energy of one energy step: phonons couple neighbouring grid points
    real(wp), parameter :: phononEnergy = Estep, gridDz = 0.25_wp
    real(wp), parameter :: eps0 = 3.0_wp, epsInf = 1.0_wp, screening = 1.0_wp
    real(wp), parameter :: cellVolume = 1.0_wp, cellLength = 1.0_wp, sigma = 0.0_wp
    real(wp), parameter :: coupling(1) = 10.0_wp
    integer, parameter :: maxScbaIterations = 10
    real(wp), parameter :: scbaTolerance = 1.0e-7_wp

    type(Tnegf) :: negf
    type(lnParams) :: params
    type(MPI_Comm) :: cartComm, kComm, enComm
    real(wp) :: coords(3,nAtoms), kPoints(3,1), kWeights(1)
    integer :: basisToMatrix(nCentral+1), ii, err
    integer, allocatable :: surfStart(:), surfEnd(:), contEnd(:), plEnd(:), cblk(:)

    surfStart = [61, 81]
    surfEnd = [60, 80]
    contEnd = [80, 100]
    plEnd = [10, 20, 30, 40, 50, 60]
    cblk = [6, 1]

    coords = 0.0_wp
    coords(3,:) = [(siteSpacing * (ii-1), ii = 1, nAtoms)]
    basisToMatrix = [(ii, ii = 1, nCentral+1)]
    kPoints = 0.0_wp
    kWeights = 1.0_wp

    call init_negf(negf)
    call set_cartesian_bare_comms(negf, comm, 1, cartComm, kComm, enComm)

    call read_HS(negf, "H_real.dat", "H_imm.dat", 0)
    call set_S_id(negf, nAtoms)

    call init_contacts(negf, 2)
    call init_structure(negf, 2, surfStart, surfEnd, contEnd, 6, plEnd, cblk)
    call init_basis(negf, coords, nCentral, basisToMatrix)
    call set_kpoints(negf, kPoints, kWeights, [1], 0)

    call get_params(negf, params)
    params%verbose = 0
    params%Emin = Emin
    params%Emax = Emin + (nPoints-1) * Estep
    params%Estep = Estep
    params%mu(1:2) = [-0.5_wp, 0.5_wp]
    params%kbT_t(1:2) = kT
    call set_params(negf, params)

    call set_elph_polaroptical(negf, coupling, phononEnergy, kT, gridDz, eps0, epsInf, screening, &
        & cellVolume, cellLength, .false., sigma, maxScbaIterations, .true.)
    call set_scba_tolerances(negf, scbaTolerance, scbaTolerance)

    if (tLayerCurrents) then
      call compute_layer_current(negf)
    else
      call compute_meir_wingreen(negf)
    end if

    ! Each task holds only the energy points it has computed
    currents = negf%curr_mat
    call MPI_Allreduce(MPI_IN_PLACE, currents, size(currents), MPI_DOUBLE_PRECISION, MPI_SUM, &
        & cartComm, err)

    call destroy_negf(negf)

  end subroutine compute_currents


  subroutine fail(message)
    character(len=*), intent(in) :: message

    integer :: err

    write(*,*) 'ERROR: ', message
    call MPI_Abort(MPI_COMM_WORLD, 1, err)

  end subroutine fail

end program test
