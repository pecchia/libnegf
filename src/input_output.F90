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


module input_output
  use ln_precision
  use ln_allocation
  use mat_def
  use sparsekit_drv
  use ln_messages
  implicit none
  private

  public :: read_H, format, is_matrix_market

  type format
     character(5) :: type ! 'PETSc','UPT','MM' (Matrix Market)
     logical :: formatted ! formatted/unformatted
     character(1) :: fmt  ! U,L,F
  end type format

  contains

!> Read a sparse matrix from file units idR (real part) and idI
!! (imaginary part). If idI <= 0 there is no imaginary file and the
!! imaginary part is taken as zero.
subroutine read_H(idR,idI, zmat,fmt)
  Integer :: idR, idI 
  Type(z_CSR) :: zmat
  Type(format) :: fmt

  !locals
  Type(r_COO) :: matr, mati
  Type(z_COO) :: mat2
  Integer :: i, count, k, nnz
  Character(10) :: tmp1, tmp2, tmp3, tmp4

  if (trim(fmt%type).eq.'MM') then
     call read_H_mm(idR, idI, zmat)
     return
  end if

  if (trim(fmt%type).eq.'PETSc') then

     read (idR,*) tmp1, tmp2, tmp3, matr%nrow, matr%ncol
     read (idR,*) tmp1, tmp4, tmp3, nnz
     read (idR,*)
     read (idR,*)
     !write(*,*) 'The number of rows (Hreal) is', matr%nrow
     !write(*,*) 'The number of columns (Hreal) is', matr%ncol
     !write(*,*) 'The number of non zero elements (Hreal) is', nnz

  elseif (trim(fmt%type).eq.'UPT') then

     if (fmt%formatted) then
        read(idR,*) matr%nrow, nnz
     else
        read(idR) matr%nrow, nnz
     endif
     !write(*,*) 'The number of rows (Hreal) is', matr%nrow
     !write(*,*) 'The number of non zero elements (Hreal) is', nnz

  end if

  select case(fmt%fmt)
  case('F')
     call create(matr,matr%nrow, matr%nrow, nnz)
     call create(mati,matr%nrow, matr%nrow, nnz)
  case('L','U')
     call create(matr,matr%nrow, matr%nrow, 2*nnz-matr%nrow)
     call create(mati,matr%nrow, matr%nrow, 2*nnz-matr%nrow) 
  end select

  if (fmt%formatted) then
     if (trim(fmt%type).eq.'PETSc') then  

        k = 0
        do i=1,nnz
           k=k+1
           read(idR,*) matr%index_j(k), matr%index_i(k), matr%nzval(k)
      
           select case(fmt%fmt)
           case('U','L') 
              if(matr%index_j(k).ne. matr%index_i(k)) then
                 k = k + 1
                 matr%index_i(k) =  matr%index_j(k-1)
                 matr%index_j(k) =  matr%index_i(k-1)
                 matr%nzval(k) = matr%nzval(k-1) 
              endif
           case('F')
           end select
        enddo

     elseif (trim(fmt%type).eq.'UPT') then
        print*,'Read Formatted UPT'
        k = 0
        do i=1,nnz
           k = k + 1           
           read(idR,*) matr%index_i(k), matr%index_j(k), matr%nzval(k)
           
           select case(fmt%fmt)
           case('U','L') 
              if(matr%index_j(k).ne. matr%index_i(k)) then
                 k = k + 1
                 matr%index_i(k) =  matr%index_j(k-1)
                 matr%index_j(k) =  matr%index_i(k-1)
                 matr%nzval(k) = matr%nzval(k-1)
              endif
           case('F')
           end select
        enddo     

     endif

  else
     if (trim(fmt%type).eq.'PETSc') then  

        k = 0
        do i=1,nnz

           k = k + 1
           read(idR) matr%index_j(k), matr%index_i(k), matr%nzval(k)
          
           select case(fmt%fmt)
           case('U','L') 
              if(matr%index_j(k).ne. matr%index_i(k)) then
                 k = k + 1
                 matr%index_i(k) =  matr%index_j(k-1)
                 matr%index_j(k) =  matr%index_i(k-1)
                 matr%nzval(k) = matr%nzval(k-1)
              endif
           case('F')
           end select
        enddo

     elseif (trim(fmt%type).eq.'UPT') then
        k = 0
        do i=1,nnz
           k = k + 1
           read(idR) matr%index_i(k), matr%index_j(k), matr%nzval(k)
           
           select case(fmt%fmt)
           case('U','L') 
              if(matr%index_j(k).ne. matr%index_i(k)) then
                 k = k + 1
                 matr%index_i(k) =  matr%index_j(k-1)
                 matr%index_j(k) =  matr%index_i(k-1)
                 matr%nzval(k) = matr%nzval(k-1)
              endif
           case('F')
           end select
        enddo     

     endif

  endif


  !write(*,*) matr%index_j, matr%index_i, matr%nzval
  ! =============================================================================
  !  Imaginary part
  ! =============================================================================
  if (idI <= 0) then
     ! No imaginary file: same pattern as the real part, zero values
     ! (k holds the number of entries read, including mirrored ones)
     mati%index_i = matr%index_i
     mati%index_j = matr%index_j
     mati%nzval = 0.0_dp
  else

  if (trim(fmt%type).eq.'PETSc') then

     read (idI,*) tmp1, tmp2, tmp3, mati%nrow, mati%ncol
     read (idI,*) tmp1, tmp4, tmp3, nnz
     read (idI,*)
     read (idI,*)
     !write(*,*) 'The number of rows (Himm) is', mati%nrow
     !write(*,*) 'The number of columns (Himm) is', mati%ncol
     !write(*,*) 'The number of non zero elements (Himm) is', nnz

  elseif (trim(fmt%type).eq.'UPT') then

     if (fmt%formatted) then
        read(idI,*) matr%nrow, nnz
     else
        read(idI) matr%nrow, nnz
     endif
     !write(*,*) 'The number of rows (Himm) is', matr%nrow
     !write(*,*) 'The number of non zero elements (Hreal) is', nnz

  end if
  
  if (mati%nrow.ne.matr%nrow) then
     call error_msg('nrow Error')
  end if
  if (matr%nnz.ne.mati%nnz) then
     call error_msg('nnz Error')
  end if
  if (matr%nnz.eq.0) then
     call error_msg('nnz Error')  
  end if

  k = 0
  do i=1,nnz

     k = k + 1
     if (fmt%formatted) then   
        if (trim(fmt%type).eq.'PETSc') then    
           read(idI,*) mati%index_j(k), mati%index_i(k), mati%nzval(k)
        elseif (trim(fmt%type).eq.'UPT') then
           read(idI,*) mati%index_i(k), mati%index_j(k), mati%nzval(k)
        endif
     else
        if (trim(fmt%type).eq.'PETSc') then
           read(idI) mati%index_j(k), mati%index_i(k), mati%nzval(k)
        elseif (trim(fmt%type).eq.'UPT') then
           read(idI) mati%index_i(k), mati%index_j(k), mati%nzval(k)
        endif
     endif

     if(mati%index_j(k).ne.matr%index_j(k)) call error_msg('Index Error')
     if(mati%index_i(k).ne.matr%index_i(k)) call error_msg('Index Error')

     select case(fmt%fmt)
     case('U','L') 
        if(mati%index_j(k).ne. mati%index_i(k)) then
           k = k + 1
           mati%index_i(k) =  mati%index_j(k-1)
           mati%index_j(k) =  mati%index_i(k-1)
           mati%nzval(k) = - mati%nzval(k-1)
        endif
     case('F')
     end select
  enddo

  end if ! idI > 0
  
  count = k
   
  call create(mat2,matr%nrow,matr%nrow,count)

  count = 0
  do i = 1, k  
       count = count + 1
       mat2%index_i(count) = matr%index_i(i)
       mat2%index_j(count) = matr%index_j(i)
       mat2%nzval(count)= cmplx(matr%nzval(i), mati%nzval(i), dp)        
  enddo
  

  call destroy(matr,mati)

  !COO-CSR conversion 
  CALL create(zmat,mat2%nrow,mat2%ncol,mat2%nnz)

  zmat%nzval=(0.d0,0.d0)

  CALL coo2csr(mat2,zmat)

  call destroy(mat2)               !deallocation Hcomplex

  !write(*,*) '(readH) matrix H read'


end subroutine read_H


! -----------------------------------------------------------------
!> True if the file on (formatted) unit id starts with a Matrix Market
!! banner. The unit is rewound.
function is_matrix_market(id) result(isMM)
  integer, intent(in) :: id
  logical :: isMM

  character(256) :: line
  integer :: ios

  isMM = .false.
  read(id, '(A)', iostat=ios) line
  if (ios == 0) then
    line = to_lower(adjustl(line))
    isMM = (line(1:14) == '%%matrixmarket')
  end if
  rewind(id)

end function is_matrix_market

! -----------------------------------------------------------------
!> Read a sparse matrix in Matrix Market coordinate format
!! ('%%MatrixMarket matrix coordinate <field> <symmetry>').
!! Entries are 'row col value' (or 'row col re im' for complex field).
!! idR : file with the real part (or the whole matrix if field=complex)
!! idI : file with the imaginary part (field real), or <= 0 if absent,
!!       in which case the imaginary part is zero.
!! symmetric/hermitian/skew-symmetric files store one triangle only:
!! the other one is reconstructed here.
subroutine read_H_mm(idR, idI, zmat)
  integer, intent(in) :: idR, idI
  type(z_CSR), intent(inout) :: zmat

  type(z_COO) :: coo
  character(32) :: field, symm, fieldI, symmI
  integer :: nrow, ncol, nnz, nrowI, ncolI, nnzI
  integer :: k, n, ii, jj, ios
  real(dp) :: re, im
  logical :: mirror

  call read_mm_header(idR, field, symm, nrow, ncol, nnz)
  if (nrow /= ncol) then
    call error_msg('read_H_mm: matrix is not square')
  end if
  mirror = (trim(symm) /= 'general')

  if (idI > 0) then
    call read_mm_header(idI, fieldI, symmI, nrowI, ncolI, nnzI)
    if (trim(field) == 'complex' .or. trim(fieldI) == 'complex') then
      call error_msg('read_H_mm: separate imaginary file requires real fields')
    end if
    if (nrowI /= nrow .or. ncolI /= ncol .or. nnzI /= nnz) then
      call error_msg('read_H_mm: real and imaginary files differ in size or nnz')
    end if
    if (trim(symmI) /= trim(symm)) then
      call error_msg('read_H_mm: real and imaginary files differ in symmetry')
    end if
  end if

  if (mirror) then
    call create(coo, nrow, ncol, 2*nnz)
  else
    call create(coo, nrow, ncol, nnz)
  end if

  im = 0.0_dp
  do k = 1, nnz
    if (trim(field) == 'complex') then
      read(idR, *, iostat=ios) ii, jj, re, im
    else
      read(idR, *, iostat=ios) ii, jj, re
    end if
    if (ios /= 0) call error_msg('read_H_mm: error reading entry', k)
    if (ii < 1 .or. ii > nrow .or. jj < 1 .or. jj > ncol) then
      call error_msg('read_H_mm: index out of range at entry', k)
    end if
    coo%index_i(k) = ii
    coo%index_j(k) = jj
    coo%nzval(k) = cmplx(re, im, dp)
  end do

  if (idI > 0) then
    do k = 1, nnz
      read(idI, *, iostat=ios) ii, jj, im
      if (ios /= 0) call error_msg('read_H_mm: error reading imaginary entry', k)
      if (ii /= coo%index_i(k) .or. jj /= coo%index_j(k)) then
        call error_msg('read_H_mm: real/imaginary index mismatch at entry', k)
      end if
      coo%nzval(k) = cmplx(real(coo%nzval(k), dp), im, dp)
    end do
  end if

  n = nnz
  if (mirror) then
    do k = 1, nnz
      if (coo%index_i(k) == coo%index_j(k)) cycle
      n = n + 1
      coo%index_i(n) = coo%index_j(k)
      coo%index_j(n) = coo%index_i(k)
      select case(trim(symm))
      case('hermitian')
        coo%nzval(n) = conjg(coo%nzval(k))
      case('skew-symmetric')
        coo%nzval(n) = -coo%nzval(k)
      case default
        coo%nzval(n) = coo%nzval(k)
      end select
    end do
  end if
  coo%nnz = n

  call create(zmat, nrow, ncol, n)
  zmat%nzval = (0.0_dp, 0.0_dp)
  call coo2csr(coo, zmat)
  call destroy(coo)

end subroutine read_H_mm

! -----------------------------------------------------------------
!> Parse the Matrix Market banner, skip comments, read the size line
subroutine read_mm_header(id, field, symm, nrow, ncol, nnz)
  integer, intent(in) :: id
  character(*), intent(out) :: field, symm
  integer, intent(out) :: nrow, ncol, nnz

  character(1024) :: line
  character(32) :: banner, object, layout
  integer :: ios

  read(id, '(A)', iostat=ios) line
  if (ios /= 0) call error_msg('read_mm_header: cannot read banner')
  read(line, *, iostat=ios) banner, object, layout, field, symm
  if (ios /= 0) call error_msg('read_mm_header: malformed banner: '//trim(line))
  banner = to_lower(banner)
  object = to_lower(object)
  layout = to_lower(layout)
  field = to_lower(field)
  symm = to_lower(symm)

  if (trim(banner) /= '%%matrixmarket' .or. trim(object) /= 'matrix') then
    call error_msg('read_mm_header: not a Matrix Market matrix file')
  end if
  if (trim(layout) /= 'coordinate') then
    call error_msg('read_mm_header: only coordinate (sparse) format is supported')
  end if
  select case(trim(field))
  case('real', 'double', 'integer', 'complex')
  case default
    call error_msg('read_mm_header: unsupported field '//trim(field))
  end select
  select case(trim(symm))
  case('general', 'symmetric', 'hermitian', 'skew-symmetric')
  case default
    call error_msg('read_mm_header: unsupported symmetry '//trim(symm))
  end select

  do
    read(id, '(A)', iostat=ios) line
    if (ios /= 0) call error_msg('read_mm_header: missing size line')
    line = adjustl(line)
    if (len_trim(line) == 0 .or. line(1:1) == '%') cycle
    exit
  end do
  read(line, *, iostat=ios) nrow, ncol, nnz
  if (ios /= 0) call error_msg('read_mm_header: malformed size line: '//trim(line))

end subroutine read_mm_header

! -----------------------------------------------------------------
pure function to_lower(str) result(low)
  character(*), intent(in) :: str
  character(len(str)) :: low
  integer :: i, ic

  low = str
  do i = 1, len(str)
    ic = iachar(str(i:i))
    if (ic >= iachar('A') .and. ic <= iachar('Z')) low(i:i) = achar(ic + 32)
  end do

end function to_lower

! -----------------------------------------------------------------




end module input_output
