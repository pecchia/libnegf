
  public :: sigma_ph_n
  public :: sigma_ph_p
  public :: sigma_ph_r
  public :: sigma_ph_r_z
  public :: complete_sigma_ph_r
  public :: check_sigma_ph_r
  public :: check_Gl_Gr
  public :: create_scratch
  public :: destroy_scratch

  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: GGn
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: GGp
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: GGr
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: Sgn
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: Sgp
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: Sgr
  logical, parameter :: memory = .true.


  subroutine create_scratch(nbl, npoints)
    integer :: nbl, npoints

    integer :: i,j,k,err

    call destroy_scratch(nbl, npoints)

    ALLOCATE(GGn(nbl,nbl,npoints),stat=err)
    ALLOCATE(GGr(nbl,nbl,npoints),stat=err)
    ALLOCATE(GGp(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgn(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgr(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgp(nbl,nbl,npoints),stat=err)

    ! Initialize everything to 0
    do k=1,npoints
      do j=1,nbl
        do i=1,nbl
          GGn(i,j,k)%nrow=0
          GGn(i,j,k)%ncol=0
          GGr(i,j,k)%nrow=0
          GGr(i,j,k)%ncol=0
          GGp(i,j,k)%nrow=0
          GGp(i,j,k)%ncol=0
          Sgn(i,j,k)%nrow=0
          Sgn(i,j,k)%ncol=0
          Sgr(i,j,k)%nrow=0
          Sgr(i,j,k)%ncol=0
          Sgp(i,j,k)%nrow=0
          Sgp(i,j,k)%ncol=0
        end do
      end do
    end do

    if (err.ne.0) then
      STOP 'ERROR: Cannot allocate GG'
    endif
    print*,'Created memory scratch',nbl,'x',nbl,'x',npoints

  end subroutine create_scratch

  subroutine destroy_scratch(nbl, npoints)
    integer :: nbl, npoints

    integer :: i, j, iE,  err

    err = 0


    do i = 1, nbl
      do j = 1, nbl
        do iE = 1, npoints

          if (allocated(GGn)) then
            if (allocated(GGn(i,j,iE)%val)) call destroy(GGn(i,j,iE))
          endif
          if (allocated(GGr)) then
            if (allocated(GGr(i,j,iE)%val)) call destroy(GGr(i,j,iE))
          endif
          if (allocated(GGp)) then
            if (allocated(GGp(i,j,iE)%val)) call destroy(GGp(i,j,iE))
          endif
          if (allocated(Sgn)) then
            if (allocated(Sgn(i,j,iE)%val)) call destroy(Sgn(i,j,iE))
          endif
          if (allocated(Sgr)) then
            if (allocated(Sgr(i,j,iE)%val)) call destroy(Sgr(i,j,iE))
          endif
          if (allocated(Sgp)) then
            if (allocated(Sgp(i,j,iE)%val)) call destroy(Sgp(i,j,iE))
          endif

        end do
      end do
    end do

    if (allocated(GGn)) DEALLOCATE(GGn,stat=err)
    if (allocated(GGr)) DEALLOCATE(GGr,stat=err)
    if (allocated(GGp)) DEALLOCATE(GGp,stat=err)
    if (allocated(Sgn)) DEALLOCATE(Sgn,stat=err)
    if (allocated(Sgr)) DEALLOCATE(Sgr,stat=err)
    if (allocated(Sgp)) DEALLOCATE(Sgp,stat=err)

    if (err.ne.0) then
      STOP 'ERROR: Cannot deallocate GG'
    endif

  end subroutine destroy_scratch

  ! ----------------------------------------------------------
  ! Search points for interpolations.
  ! Wq has a sign (+/- Wq)
  !
  subroutine search_points(negf, Wq, Epnt, i1, i2, E1, E2)
    use energy_mesh, only : elem
    type(TNegf) :: negf
    real(dp) :: Wq
    real(dp), dimension(:), allocatable :: Epnt
    integer, intent(out) :: i1,i2
    real(dp), intent(out) :: E1, E2

    integer :: iE, iel, istart, iend, ip
    Type(elem), pointer :: pel
    real(dp) :: En

    if (Wq.eq.0) then
      i1 = negf%iE
      i2 = negf%iE
      E1 = real(negf%Epnt)
      E2 = real(negf%Epnt)
      return
    endif

    if (allocated(Epnt)) then
      ! Remove offset such that search can work on Epnt(1..N)
      iE = negf%iE - negf%Np_n(1) - negf%Np_n(2) - negf%n_poles
      En = real(negf%Epnt) + Wq !Wq carry the right sign
      !print*
      !print*,'iE', iE, real(negf%Epnt) + Wq

      if (sign(1.0_dp,Wq) .gt. 0) then
        i2 = iE + 1
        if (i2.gt.size(Epnt)) then
          i2 = size(Epnt)
        else
          do while (Epnt(i2) .lt. En)
            i2 = i2 + 1
          end do
        endif
        i1 = i2 - 1
      else
        i1 = iE - 1
        if (i1.lt.1) then
          i1 = 1
        else
          do while (Epnt(i1) .gt. En)
            i1 = i1 - 1
          end do
        endif
        i2 = i1 + 1
      endif
      E1 = Epnt(i1)
      E2 = Epnt(i2)
      ! add back the offset to the point
      i1 = i1 + negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
      i2 = i2 + negf%Np_n(1) + negf%Np_n(2) + negf%n_poles

    else

      !if (.not.allocated(negf%emesh%pactive)) STOP 'emesh not initialized'

      En = real(negf%Epnt) + Wq

      if (sign(1.0_dp,Wq) .gt. 0) then
        istart = negf%emesh%iactive
        iend = negf%emesh%maxind

        elloop1: do iel = istart, iend
          pel => negf%emesh%pactive(iel)%pelem
          do ip = 1, 3
            if (pel%pnt(ip) .gt. En) then
              exit elloop1
            endif
          end do
        end do elloop1
        i1 = pel%map(ip-1)
        i2 = pel%map(ip)
        E1 = pel%pnt(ip-1)
        E2 = pel%pnt(ip)
      else
        istart = negf%emesh%iactive
        iend = 1

        elloop2: do iel = istart, iend, -1
          pel => negf%emesh%pactive(iel)%pelem
          do ip = 3, 1, -1
            if (pel%pnt(ip) .lt. En) then
              exit elloop2
            endif
          end do
        end do elloop2
        i1 = pel%map(ip)
        i2 = pel%map(ip+1)
        E1 = pel%pnt(ip)
        E2 = pel%pnt(ip+1)
      end if
    end if
    !print*
    !print*,E1,En,E2
    !print*,'interpolate between:',i1,i2

  end subroutine search_points


  subroutine interpolation(i1,i2, E1, E2, E, path, name, G_interp)
    integer, intent(in) :: i1, i2
    real(dp) :: E1, E2, E
    type(z_DNS), dimension(:,:) :: G_interp
    CHARACTER(*) :: path
    CHARACTER(*) :: name

    !local variables
    type(z_DNS) :: work1,work2
    integer :: i

    do i = 1, size(G_interp,1)

      call create(work1,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
      work1%val = (0.0_dp,0.0_dp)
      call read_blkmat(work1, path, name, i, i, i1)

      if (E1.ne.E2 .and. i1.ne.i2) then

        call create(work2,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
        work2%val = (0.0_dp,0.0_dp)
        call read_blkmat(work2, path, name, i, i, i2)

        G_interp(i,i)%val = ((E-E1)*work2%val + (E2-E)*work1%val)/(E2-E1)

        call destroy(work2)

      else

        G_interp(i,i)%val = work1%val

      endif

      call destroy(work1)

    end do


  end subroutine interpolation

  !--------------------------------------------------------------------------------
  subroutine Sigma_ph_r(negf,Epnt)

    type(Tnegf) :: negf
    real(dp), dimension(:),allocatable :: Epnt

    !Local variables
    type(z_DNS) :: work1,work2
    type(z_DNS) :: Mq_mat
    type(z_DNS), dimension(:,:), allocatable :: Sigma_r, G_r_interP, G_r_interN
    type(z_DNS), dimension(:,:), allocatable :: G_n_interP, G_n_interN

    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Nq
    real(dp), dimension(:), pointer :: Mq
    logical, dimension(:), pointer :: selmodes
    integer :: i, m, iE, i1, i2
    integer :: nummodes, numselmodes, nbl
    real(dp) :: E1, E2, En

    nbl = negf%str%num_PLs


    selmodes => negf%elph%selmodes
    Mq => negf%elph%Mq
    Wq => negf%elph%Wq
    Nq => negf%elph%Nq
    nummodes = negf%elph%nummodes
    numselmodes = negf%elph%numselmodes

    iE=negf%iE

    allocate(G_r_interP(nbl,nbl))
    allocate(G_r_interN(nbl,nbl))
    allocate(G_n_interP(nbl,nbl))
    allocate(G_n_interN(nbl,nbl))
    allocate(Sigma_r(nbl,nbl))

    associate(indblk=>negf%str%mat_PL_start)
    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_r(i,i), m, m)
      Sigma_r(i,i)%val = (0.0_dp, 0.0_dp)
      call create(G_r_interP(i,i), m, m)
      call create(G_r_interN(i,i), m, m)
      call create(G_n_interP(i,i), m, m)
      call create(G_n_interN(i,i), m, m)
      !call create(G_r(i,i), m, m)
      !call read_blkmat(G_r(i,i),negf%scratch_path,'G_r_',i,i,iE)
    end do
    end associate
    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,Wq(m),Epnt,i1,i2,E1,E2)
      En = real(negf%Epnt)+Wq(m)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_n_', G_n_interP)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_r_', G_r_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,-Wq(m),Epnt,i1,i2,E1,E2)
      En = real(negf%Epnt)-Wq(m)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_n_', G_n_interN)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_r_', G_r_interN)

      do i = 1, nbl

        i1 = Sigma_r(i,i)%nrow

        call create(work1,i1,i1)

        work1%val = (0.0_dp,0.0_dp)

        if (negf%elph%selfene_gr) then
          !print*,'SelfEneR_Gr'
          ! Via I: should be exact for E >> Ef_max
          work1%val = (Nq(m)+1.0_dp)*G_r_interN(i,i)%val + Nq(m)*G_r_interP(i,i)%val
          ! Via II: should be exact for E << Ef_min
          !work1%val = (Nq(m)+1.0_dp)*G_r_interP(i,i)%val + Nq(m)*G_r_interN(i,i)%val
          ! Via III: should work as a compromise
          !work1%val = Nq(m)*(G_r_interP(i,i)%val + G_r_interN(i,i)%val) + G_r(i,i)%val
        endif

        if (negf%elph%selfene_gless) then
          !print*,'SelfEneR_G<'
          ! should be 1/2 [G<- - G<+] == i/2 [Gn- - Gn+]
          work1%val = work1%val + &
              (0.0_dp,0.5_dp) * (G_n_interN(i,i)%val - G_n_interP(i,i)%val)
        endif

        if (negf%elph%selfene_gless .or. negf%elph%selfene_gr) then

          call create_id(Mq_mat,i1,Mq(m))

          call prealloc_mult(Mq_mat, work1, work2)

          call destroy(work1)

          call prealloc_mult(work2, Mq_mat, work1)

          call destroy(work2)

          Sigma_r(i,i)%val = Sigma_r(i,i)%val + work1%val

          call destroy(work1, Mq_mat)

        else

          Sigma_r(i,i)%val = (0.0_dp,0.0_dp)
          call destroy(work1)

        endif

      end do

    end do

    ! throw away all non diagonal parts
    !do i = 1, nbl
    !    i1 = Sigma_r(i,i)%nrow
    !    do m = 1, i1
    !       Sigma_r(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp)
    !       Sigma_r(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp)
    !    end do
    ! end do
    !

    do i = 1, nbl
      !print*
      !print*,'(sigma_r) Sigma_ph_r',maxval(abs(Sigma_r(i,i)%val)), iE
      call write_blkmat(Sigma_r(i,i),negf%scratch_path,'Sigma_ph_r_',i,i,iE)
      call destroy(Sigma_r(i,i))
      call destroy(G_r_interP(i,i))
      call destroy(G_r_interN(i,i))
      call destroy(G_n_interP(i,i))
      call destroy(G_n_interN(i,i))
      !call destroy(G_r(i,i))
    end do

    deallocate(Sigma_r,G_r_interP,G_r_interN,G_n_interP,G_n_interN)
    !deallocate(G_r)

  end subroutine Sigma_ph_r


  subroutine Sigma_ph_r_z(negf,z)

    type(Tnegf) :: negf
    complex(dp), intent(in) :: z

    !Local variables
    type(z_DNS) :: work1,work2
    type(z_DNS) :: Mq_mat
    type(z_DNS), dimension(:,:), allocatable :: Sigma_r, G_r
    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Nq
    real(dp), dimension(:), pointer :: Mq
    logical, dimension(:), pointer :: selmodes
    integer :: i, m, iE,i1
    integer :: nummodes, numselmodes, nbl

    nbl = negf%str%num_PLs

    selmodes => negf%elph%selmodes
    Mq => negf%elph%Mq
    Wq => negf%elph%Wq
    Nq => negf%elph%Nq
    nummodes = negf%elph%nummodes
    numselmodes = negf%elph%numselmodes

    iE=negf%iE

    allocate(Sigma_r(nbl,nbl))
    allocate(G_r(nbl,nbl))
    
    associate(indblk => negf%str%mat_PL_start)
    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_r(i,i), m, m)
      Sigma_r(i,i)%val = (0.0_dp, 0.0_dp)
      call create(G_r(i,i), m, m)
    end do
    end associate

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle

      do i = 1, nbl

        call read_blkmat(G_r(i,i),negf%scratch_path,'G_r_',i,i,iE)

        i1 = Sigma_r(i,i)%nrow

        if (negf%elph%selfene_gr) then

          call create(work1,i1,i1)

          work1%val = (0.0_dp,0.0_dp)

          work1%val = (2.0_dp*Nq(m)+1.0_dp)*G_r(i,i)%val

          call create_id(Mq_mat,i1,Mq(m))

          call prealloc_mult(Mq_mat, work1, work2)

          call destroy(work1)

          call prealloc_mult(work2, Mq_mat, work1)

          call destroy(work2)

          Sigma_r(i,i)%val = Sigma_r(i,i)%val +  work1%val

          call destroy(work1, Mq_mat)

        else
          Sigma_r(i,i)%val = (0.0_dp,0.0_dp)
        endif

      end do

    end do

    do i = 1, nbl
      call write_blkmat(Sigma_r(i,i),negf%scratch_path,'Sigma_ph_r_',i,i,iE)
      call destroy(Sigma_r(i,i))
      call destroy(G_r(i,i))
    end do

  end subroutine Sigma_ph_r_z
 

  !*********************************************************************
  !ADD Hilbert-transform part to Sigma_ph_r
  !Need to set back FFT transforms
  !********************************************************************
  subroutine complete_sigma_ph_r(negf, Epnt, ioffset)

    type(Tnegf) :: negf
    real(dp), dimension(:) :: Epnt
    integer :: ioffset

    ! Locals
    integer :: i,j,n,m, iE, bl, sizebl,i_start,i_stop, ierr, nummodes
    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Mq
    real(dp) :: Wmax, dE, tmp

    character(4) :: ofblki
    complex(dp), dimension(:), allocatable :: temp1
    type(z_DNS), dimension(:), allocatable :: Gn_E, Sigma_r_E
    !type(z_DNS) :: work1, work2, work3, Mq_mat
    logical, dimension(:), pointer :: selmodes

    Mq => negf%elph%Mq
    Wq => negf%elph%Wq
    selmodes => negf%elph%selmodes

    n = size(Epnt)
    nummodes = size(Wq)

    Wmax = maxval(Wq)
    dE = Epnt(2)-Epnt(1)
    ! ENERGY-INTERVAL FOR SELF-ENERGIES:
    i_start = int(aint(Wmax/dE) + 1)
    i_stop =  (n - i_start)
    i_start = i_start + 1

    ! PROBLEM:
    ! The Hilb-Transf should be resticted on the sub interval
    ! Emin+m*Wmax ... Emax-m*Wmax
    ! However the number of points should be always 2^p
    ! This condition is difficult to fulfill.
    ! Possible solution: interpolation between grids.

    ! CREATE THE ARRAY OF MATRICES (one for each E-point)
    allocate(Gn_E(n),stat=ierr)
    allocate(Sigma_r_E(n),stat=ierr)
    if (ierr.ne.0) STOP 'ERROR in allocation of Gn_E'
    call log_allocate(temp1,n)

    print*
    print*,'HILBERT TRANSFORM memory:',sizebl*sizebl*n*16
    print*,'HILBERT TRANSFORM interval:',i_start,i_stop

    ! LOOP ON blocks
    do bl = 1, negf%str%num_PLs

      sizebl =  negf%str%mat_PL_start(bl+1) - negf%str%mat_PL_start(bl)
      if (bl.le.9999) write(ofblki,'(i4.4)') bl
      if (bl.gt.9999) stop 'ERROR: too many blks (> 9999)'

      ! LOAD ALL G< FROM FILES and store in Gn_E(iE)
      do iE = 1, n
        call create(Gn_E(iE),sizebl,sizebl)
        call read_blkmat(Gn_E(iE), negf%scratch_path, 'G_n_', bl, bl,iE+ioffset)
      end do

      ! LOAD ALL Sigma_r in the right interval
      do iE = i_start, i_stop
        call create(Sigma_r_E(iE),sizebl,sizebl)
        call read_blkmat(Sigma_r_E(iE), negf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
      end do

      do m = 1, nummodes
        if (.not.selmodes(m)) cycle

        ! COMPUTE   Mq G< Mq   (assume diagonal now)
        do iE = 1, n
          !call create_id(Mq_mat,sizebl,Mq(m))
          !call prealloc_mult(Mq_mat, Gn_E(iE), work1)
          !call prealloc_mult(work1, Mq_mat, work2)
          !call destroy(work1, Mq_mat)
          !Gn_E(iE)%val = work2%val
          !call destroy(work2)
          tmp = Mq(m)*Mq(m)
          Gn_E(iE)%val = tmp*Gn_E(iE)%val
        end do

        ! PERFORM Hilbert in Energy for each (i,j)
        do j = 1, sizebl
          do i = 1, sizebl

            !          k = (j-1)*sizebl+i
            !          do iE = 1, n
            !             if (iE.le.99999) write(ofpnt,'(i5.5)') iE+ioffset
            !             filename = 'G_n_'//trim(ofblki)//'_'//trim(ofblki)//'_'//trim(ofpnt)//'.dat'
            !             open(100,file=trim(negf%scratch_path)//trim(filename), access='DIRECT', recl=4)
            !             READ(100, rec = k)  temp1(iE)
            !             close(100)
            !          end do

            ! SETUP a vector out of all G<_ij(E)
            ! Here we could perform an efficient interpolation on a regular grid 2^p
            do iE = 1, n
              temp1(iE) = Gn_E(iE)%val(i,j)
            end do

            Wq(m) = Wq(m)*2.0_dp*pi/((n-1)*dE)

            !call Hilbert_shift(temp1, Wq(m))

            Wq(m) = Wq(m)*((n-1)*dE)/(2.0_dp*pi)

            ! UPDATE the self-energies with the Hilbert part.
            do iE = i_start, i_stop
              ! Should be  -i/2 H[ G<+ - G<-  ] = 1/2 H[ Gn+ - Gn- ]
              Sigma_r_E(iE)%val(i,j) =  Sigma_r_E(iE)%val(i,j) + (0.5_dp, 0.0)* temp1(iE)
            end do


            !          do iE = i_start+1, i_stop
            !             if (iE.le.99999) write(ofpnt,'(i5.5)') iE+ioffset
            !             filename = 'Sigma_ph_r_'//trim(ofblki)//'_'//trim(ofblki)//'_'//trim(ofpnt)//'.dat'
            !             open(100,file=trim(negf%scratch_path)//trim(filename), access='DIRECT', recl=4)
            !             READ (100, rec = k) temp2
            !             temp2 = temp2 - (0.0_dp, 0.5_dp)* temp1(iE)
            !             WRITE (100, rec = k) temp2
            !             close(100)
            !          end do

          end do ! Loop on block size
        end do ! Loop on block size

      end do !Loop on modes

      do iE = i_start, i_stop
        call write_blkmat(Sigma_r_E(iE), negf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
      end do

      do iE = 1, n
        call destroy(Gn_E(iE))
      end do
      do iE = i_start, i_stop
        call destroy(Sigma_r_E(iE))
      end do

    end do !Loop on blocks

    call log_deallocate(temp1)
    deallocate(Gn_E)
    deallocate(Sigma_r_E)

  end subroutine complete_sigma_ph_r

  ! ******************************************************************************
  ! Computes Sigma_ph_n and save it file
  ! ******************************************************************************
  subroutine Sigma_ph_n(negf,Epnt)

    type(Tnegf) :: negf
    real(dp), dimension(:),allocatable :: Epnt

    !Local variables
    type(z_DNS) :: work1,work2
    type(z_DNS) :: Mq_mat
    type(z_DNS), dimension(:,:), allocatable :: Sigma_n, G_n_interP, G_n_interN

    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Nq
    real(dp), dimension(:), pointer :: Mq
    logical, dimension(:), pointer :: selmodes

    integer :: i, m, iE, i1,i2
    integer :: nummodes, numselmodes, nbl
    real(dp) :: E1, E2, En

    nbl = negf%str%num_PLs

    selmodes => negf%elph%selmodes
    Wq => negf%elph%Wq
    Mq => negf%elph%Mq
    Nq => negf%elph%Nq
    numselmodes = negf%elph%numselmodes
    nummodes = negf%elph%nummodes

    iE = negf%iE

    allocate(G_n_interP(nbl,nbl))
    allocate(G_n_interN(nbl,nbl))
    allocate(Sigma_n(nbl,nbl))

    associate(indblk=>negf%str%mat_PL_start)
    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_n(i,i), m, m)
      Sigma_n(i,i)%val = (0.0_dp, 0.0_dp)
      call create(G_n_interP(i,i), m, m)
      call create(G_n_interN(i,i), m, m)
    end do
    end associate

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle
      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,Wq(m),Epnt,i1,i2,E1,E2)

      En = real(negf%Epnt)+Wq(m)

      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_n_', G_n_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,-Wq(m),Epnt,i1,i2,E1,E2)

      En = real(negf%Epnt)-Wq(m)

      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_n_', G_n_interN)


      do i = 1, nbl

        i1 =  G_n_interN(i,i)%nrow

        call create(work1, i1, i1)

        work1%val = (Nq(m)+1.0_dp)*G_n_interP(i,i)%val + Nq(m)* G_n_interN(i,i)%val

        call create_id(Mq_mat,i1,Mq(m))

        call prealloc_mult( Mq_mat, work1, work2)

        call destroy(work1)

        call prealloc_mult( work2, Mq_mat, work1)

        call destroy(work2)

        Sigma_n(i,i)%val = Sigma_n(i,i)%val + work1%val

        call destroy(work1, Mq_mat)

      end do

    end do

    !Throw away all non diagonal parts
    !do i = 1, nbl
    !   i1 = Sigma_n(i,i)%nrow
    !   do m = 1, i1
    !      Sigma_n(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp)
    !      Sigma_n(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp)
    !   end do
    !end do
    !

    do i = 1, nbl
      call write_blkmat(Sigma_n(i,i),negf%scratch_path,'Sigma_ph_n_',i,i,iE)
      call destroy(Sigma_n(i,i))
      call destroy(G_n_interN(i,i))
      call destroy(G_n_interP(i,i))
    end do

    deallocate(Sigma_n,G_n_interN, G_n_interP)

  end subroutine Sigma_ph_n


  ! ******************************************************************************
  ! Computes Sigma_ph> and save it file
  ! ******************************************************************************
  subroutine Sigma_ph_p(negf,Epnt)

    type(Tnegf) :: negf
    real(dp), dimension(:),allocatable :: Epnt

    !Local variables
    type(z_DNS) :: work1,work2
    type(z_DNS) :: Mq_mat
    type(z_DNS), dimension(:,:), allocatable :: Sigma_p, G_p_interP, G_p_interN

    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Nq
    real(dp), dimension(:), pointer :: Mq
    logical, dimension(:), pointer :: selmodes

    integer :: i, m, iE, i1,i2
    integer :: nummodes, numselmodes, nbl
    real(dp) :: E1, E2, En

    nbl = negf%str%num_PLs

    selmodes => negf%elph%selmodes
    Wq => negf%elph%Wq
    Mq => negf%elph%Mq
    Nq => negf%elph%Nq
    numselmodes = negf%elph%numselmodes
    nummodes = negf%elph%nummodes

    iE = negf%iE

    allocate(G_p_interP(nbl,nbl))
    allocate(G_p_interN(nbl,nbl))
    allocate(Sigma_p(nbl,nbl))

    associate(indblk=>negf%str%mat_PL_start)
    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_p(i,i), m, m)
      Sigma_p(i,i)%val = (0.0_dp, 0.0_dp)
      call create(G_p_interP(i,i), m, m)
      call create(G_p_interN(i,i), m, m)
    end do
    end associate

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle
      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,Wq(m),Epnt,i1,i2,E1,E2)

      En = real(negf%Epnt)+Wq(m)

      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_p_', G_p_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,-Wq(m),Epnt,i1,i2,E1,E2)
      En = real(negf%Epnt)-Wq(m)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_p_', G_p_interN)

      do i = 1, nbl

        i1 =  G_p_interN(i,i)%nrow

        call create(work1, i1, i1)

        work1%val = (Nq(m)+1.0_dp)*G_p_interN(i,i)%val + Nq(m)* G_p_interP(i,i)%val

        call create_id(Mq_mat,i1,Mq(m))

        call prealloc_mult( Mq_mat, work1, work2)

        call destroy(work1)

        call prealloc_mult ( work2, Mq_mat, work1)

        call destroy(work2)

        Sigma_p(i,i)%val = Sigma_p(i,i)%val + work1%val

        call destroy(work1, Mq_mat)

      end do

    end do

    ! throw away all non diagonal parts
    !do i = 1, nbl
    !   i1 = Sigma_p(i,i)%nrow
    !   do m = 1, i1
    !      Sigma_p(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp)
    !      Sigma_p(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp)
    !   end do
    !end do
    !

    do i = 1, nbl
      call write_blkmat(Sigma_p(i,i),negf%scratch_path,'Sigma_ph_p_',i,i,iE)
      call destroy(Sigma_p(i,i))
      call destroy(G_p_interN(i,i))
      call destroy(G_p_interP(i,i))
    end do

    deallocate(Sigma_p,G_p_interN, G_p_interP)

  end subroutine Sigma_ph_p
  
  ! ----------------------------------------------------------
  subroutine check_Gl_Gr(negf)
    type(Tnegf) :: negf
    real(dp), dimension(:),allocatable :: Epnt
    integer :: ioffset

    type(z_DNS), dimension(:,:), allocatable :: G_r, G_p, G_n
    type(z_DNS) :: A, T
    integer :: nbl, n, sizebl, i_start, i_stop, psize, maxpos(2)
    real(dp) :: Wmax, maxdev, tmp, maxG

    nbl = negf%str%num_PLs

    allocate(G_n(nbl,nbl))
    allocate(G_p(nbl,nbl))
    allocate(G_r(nbl,nbl))

    maxdev = 0.0_dp
    maxG = 0.0_dp
    psize = 0

    do n = 1, nbl
      sizebl = negf%str%mat_PL_start(n+1)-negf%str%mat_PL_start(n)
      call create(G_r(n,n),sizebl,sizebl)
      call create(G_n(n,n),sizebl,sizebl)
      call create(G_p(n,n),sizebl,sizebl)
      call read_blkmat(G_r(n,n), negf%scratch_path, 'G_r_',n,n, negf%iE)
      call read_blkmat(G_n(n,n),negf%scratch_path,'G_n_',n,n,negf%iE)
      call read_blkmat(G_p(n,n),negf%scratch_path,'G_p_',n,n,negf%iE)

      call zspectral(G_r(n,n),G_r(n,n),0,A)

      call create(T,sizebl,sizebl)


      T%val = G_n(n,n)%val + G_p(n,n)%val - A%val

      tmp = maxval(abs(A%val))
      if (tmp .gt. maxG) maxG=tmp

      tmp = maxval(abs(T%val))/maxval(abs(A%val))
      if (tmp .gt. maxdev) then
        maxdev = tmp
        maxpos = maxloc(abs(T%val)) + psize
      endif

      psize = psize + sizebl

      call destroy(G_r(n,n))
      call destroy(G_n(n,n))
      call destroy(G_p(n,n))
      call destroy(A)
      call destroy(T)

    end do

    !print*,'CHECK Gn+Gp=Gr-Ga',negf%iE, maxG, maxdev

    deallocate(G_n, G_p, G_r)

  end subroutine check_Gl_Gr

  ! ----------------------------------------------------------
  subroutine check_sigma_ph_r(negf)
    type(Tnegf) :: negf

    type(z_DNS), dimension(:,:), allocatable :: Sigma_r, Sigma_p, Sigma_n
    type(z_DNS) :: Gam, T
    integer :: nbl, n, sizebl, psize, maxpos(2)
    real(dp) :: maxdev, tmp, maxG

    nbl = negf%str%num_PLs

    allocate(Sigma_n(nbl,nbl))
    allocate(Sigma_p(nbl,nbl))
    allocate(Sigma_r(nbl,nbl))

    maxdev = 0.0_dp
    maxG = 0.0_dp
    psize = 0

    do n = 1, nbl
      sizebl = negf%str%mat_PL_start(n+1)-negf%str%mat_PL_start(n)
      call create(Sigma_r(n,n),sizebl,sizebl)
      call create(Sigma_n(n,n),sizebl,sizebl)
      call create(Sigma_p(n,n),sizebl,sizebl)
      call read_blkmat(Sigma_r(n,n), negf%scratch_path, 'Sigma_ph_r_',n,n, negf%iE)
      call read_blkmat(Sigma_n(n,n),negf%scratch_path,'Sigma_ph_n_',n,n,negf%iE)
      call read_blkmat(Sigma_p(n,n),negf%scratch_path,'Sigma_ph_p_',n,n,negf%iE)

      call zspectral(Sigma_r(n,n),Sigma_r(n,n),0,Gam)

      call create(T,sizebl,sizebl)

      !print*,'CHECK Sigma_ph_n+Sigma_ph_p',maxval(abs(Sigma_n(n,n)%val+Sigma_p(n,n)%val))

      T%val = Sigma_n(n,n)%val + Sigma_p(n,n)%val - Gam%val

      tmp = maxval(abs(Gam%val))
      if (tmp .gt. maxG) maxG=tmp

      tmp = maxval(abs(T%val))/maxval(abs(Gam%val))
      if (tmp .gt. maxdev) then
        maxdev = tmp
        maxpos = maxloc(abs(T%val)) + psize
      endif

      psize = psize + sizebl

      call destroy(Sigma_r(n,n))
      call destroy(Sigma_n(n,n))
      call destroy(Sigma_p(n,n))
      call destroy(Gam)
      call destroy(T)

    end do

    !print*,'CHECK Sigma_ph_r',negf%iE, maxG, maxdev

    deallocate(Sigma_n, Sigma_p, Sigma_r)


  end subroutine check_sigma_ph_r


  ! READ Matrices
  subroutine read_blkmat(Matrix, path, name, i, j, iE)
    type(z_DNS), intent(inout) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    integer, intent(in) :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    logical :: lex
    integer :: i1,i2
    complex(dp) :: mat_el

    if (memory) then
      select case(trim(name))
      case('G_n_')
        Matrix%val = GGn(i,j,iE)%val
      case('G_r_')
        Matrix%val = GGr(i,j,iE)%val
      case('G_p_')
        Matrix%val = GGp(i,j,iE)%val
      case('Sigma_ph_r_')
        Matrix%val = Sgr(i,j,iE)%val
      case('Sigma_ph_n_')
        Matrix%val = Sgn(i,j,iE)%val
      case('Sigma_ph_p_')
        Matrix%val = Sgp(i,j,iE)%val
      case default
        stop 'internal error: read_blkmat does not correspond'
      end select
      return
    endif

    Matrix%val = (0.0_dp,0.0_dp)

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'

    inquire(file=trim(path)//trim(filename),EXIST=lex)
    if (.not.lex) then
      RETURN
      !WRITE(*,*) 'ERROR: FILE '//trim(filename)//' doES NOT EXIST'
      !STOP
    endif

    open(9091,file=trim(path)//trim(filename), access='STREAM')

    call inmat_c(9091,.false.,Matrix%val,Matrix%nrow,Matrix%ncol)

    !open(9001,file=trim(path)//trim(filename), access='DIRECT', recl=4)
    !call direct_in_c(9001,Matrix%val,Matrix%nrow)

    close(9091)

  end subroutine read_blkmat

  ! WRITE Matrices
  subroutine write_blkmat(Matrix, path, name, i, j, iE)
    type(z_DNS), intent(in) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    integer, intent(in) :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    integer :: m,n

    if (memory) then

      select case(trim(name))
      case('G_n_')
        if (.not.allocated(GGn(i,j,iE)%val)) then
          call create(GGn(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGn(i,j,iE)%val = Matrix%val
      case('G_r_')
        if (.not.allocated(GGr(i,j,iE)%val)) then
          call create(GGr(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGr(i,j,iE)%val = Matrix%val
      case('G_p_')
        if (.not.allocated(GGp(i,j,iE)%val)) then
          call create(GGp(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGp(i,j,iE)%val = Matrix%val
      case('Sigma_ph_r_')
        if (.not.allocated(Sgr(i,j,iE)%val)) then
          call create(Sgr(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgr(i,j,iE)%val = Matrix%val
      case('Sigma_ph_n_')
        if (.not.allocated(Sgn(i,j,iE)%val)) then
          call create(Sgn(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgn(i,j,iE)%val = Matrix%val
      case('Sigma_ph_p_')
        if (.not.allocated(Sgp(i,j,iE)%val)) then
          call create(Sgp(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgp(i,j,iE)%val = Matrix%val
      case default
        stop 'internal error: write_blkmat does not correspond'
      end select
      return
    endif

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'

    open(9001,file=trim(path)//trim(filename), access='STREAM', status='REPLACE')

    call outmat_c(9001,.false.,Matrix%val,Matrix%nrow,Matrix%ncol) !,1.0d-36)

    !open(9001,file=trim(path)//trim(filename), status='REPLACE', access='DIRECT', recl=4)
    !call direct_out_c(9001,Matrix%val,Matrix%nrow)

    close(9001)

  end subroutine write_blkmat
  !****************************************************************************
  !
  ! Calculate G_p=iG> contributions due to el-ph
  ! Writing on memory
  !
  !****************************************************************************
!!$  subroutine calculate_Gp_ph(negf,ESH,iter,Gp)
!!$
!!$    type(Tnegf) :: negf
!!$    type(z_DNS), dimension(:,:) :: ESH, Gp
!!$    integer :: iter
!!$
!!$    Type(z_DNS), dimension(:,:), allocatable :: Sigma_ph_p
!!$    Type(z_DNS) :: Ga, work1, work2
!!$    integer :: n, k, nbl, nrow, ierr
!!$
!!$    nbl = negf%str%num_PLs
!!$    ALLOCATE(Sigma_ph_p(nbl,nbl),stat=ierr)
!!$    if (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'
!!$
!!$
!!$    do n = 1, nbl
!!$
!!$      nrow = ESH(n,n)%nrow
!!$
!!$      call create(Sigma_ph_p(n,n), nrow, nrow)
!!$
!!$      Sigma_ph_p(n,n)%val = (0.0_dp, 0.0_dp)
!!$      if (iter .gt. 0) then
!!$         call read_blkmat(Sigma_ph_p(n,n),negf%scratch_path,'Sigma_ph_p_',n,n,negf%iE)
!!$      else
!!$         call write_blkmat(Sigma_ph_p(n,n),negf%scratch_path,'Sigma_ph_p_',n,n,negf%iE)
!!$      endif
!!$
!!$    end do
!!$
!!$    do n = 1, nbl
!!$
!!$      do k = 1, nbl
!!$
!!$         if (Gr(n,k)%nrow.gt.0) then
!!$            call zdagger(Gr(n,k),Ga)
!!$            call prealloc_mult(Gr(n,k), Sigma_ph_p(k,k), work1)
!!$            call prealloc_mult(work1, Ga, work2)
!!$            Gp(n,n)%val = Gp(n,n)%val + work2%val
!!$            call destroy(work1, work2, Ga)
!!$         endif
!!$
!!$      end do
!!$
!!$    end do
!!$
!!$    do n = 1, nbl
!!$      call destroy(Sigma_ph_p(n,n))
!!$    end do
!!$
!!$    DEALLOCATE(Sigma_ph_p)
!!$
!!$  end subroutine calculate_Gp_ph
