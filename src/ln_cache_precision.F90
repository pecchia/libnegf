!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute and/or modify it          !
!! under the terms of the GNU Lesser General Public License as published    !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!

!> Precision of the matrix caches (G_r, G_n, SGF, Sigma of the inelastic
!> interactions) and of the self-energy kernels working on them.
!> Selected at compile time: fypp -DSP_CACHE (CMake WITH_SP_CACHE) gives
!> single precision, the default is double precision.
!>   x_DNS          : dense block type stored in the caches
!>   lp             : real/complex kind of the kernels
!>   MPI_XX_COMPLEX : MPI datatype matching complex(lp)
module ln_cache_precision
  use ln_precision, only : sp, dp
#:if defined("SP_CACHE")
  use mat_def, only : x_DNS => c_DNS
  use mpi_f08, only : MPI_XX_COMPLEX => MPI_COMPLEX
#:else
  use mat_def, only : x_DNS => z_DNS
  use mpi_f08, only : MPI_XX_COMPLEX => MPI_DOUBLE_COMPLEX
#:endif
  implicit none
  private

  public :: x_DNS, lp, MPI_XX_COMPLEX

#:if defined("SP_CACHE")
  integer, parameter :: lp = sp
#:else
  integer, parameter :: lp = dp
#:endif

end module ln_cache_precision
