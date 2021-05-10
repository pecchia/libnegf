module cudautils
   use ln_precision
   use openacc
   use cublas_v2
   use cusolverDn
   use mat_def
   implicit none
   private
   public :: createGPU
   public :: copyToGPU
   public :: copyFromGPU
   public :: deleteGPU
   public :: createAll
   public :: destroyAll

   public :: copyMatToGPU
   public :: copyMatFromGPU
   public :: deleteMatGPU
   public :: createMatGPU
   public :: copy_ESH_toGPU
   
   public :: copy_mat_gpu
   public :: matmul_gpu
   public :: matmul_gpu_dag
   public :: add_gpu
   public :: add_cublas
   public :: sum_gpu
   public :: inverse_gpu
   public :: init_mat_gpu
   public :: zspectral_gpu

   interface copyToGPU
      module procedure copyToGPU_sp
      module procedure copyToGPU_dp
   end interface copyToGPU

   interface createGPU
      module procedure createGPU_sp
      module procedure createGPU_dp
   end interface createGPU

   interface copyFromGPU
      module procedure copyFromGPU_sp
      module procedure copyFromGPU_dp
   end interface copyFromGPU

   interface deleteGPU
      module procedure deleteGPU_sp
      module procedure deleteGPU_dp
   end interface deleteGPU

   interface createAll
      module procedure createAll_sp
      module procedure createAll_dp
   end interface createAll

   interface destroyAll
      module procedure destroyAll_sp
      module procedure destroyAll_dp
   end interface destroyAll


   
   interface copyMatToGPU
      module procedure copyMatToGPU_sp
      module procedure copyMatToGPU_dp
   end interface copyMatToGPU

   interface copyMatFromGPU
      module procedure copyMatFromGPU_sp
      module procedure copyMatFromGPU_dp
   end interface copyMatFromGPU

   interface createMatGPU
      module procedure createMatGPU_sp
      module procedure createMatGPU_dp
   end interface createMatGPU

   interface deleteMatGPU
      module procedure deleteMatGPU_sp
      module procedure deleteMatGPU_dp
   end interface deleteMatGPU

   interface copy_ESH_toGPU
      module procedure copy_ESH_toGPU_sp
      module procedure copy_ESH_toGPU_dp
   end interface copy_ESH_toGPU


   
   interface copy_mat_gpu
      module procedure copy_mat_gpu_sp
      module procedure copy_mat_gpu_dp
   end interface copy_mat_gpu

   interface matmul_gpu
      module procedure matmul_gpu_sp
      module procedure matmul_gpu_dp
   end interface matmul_gpu

   interface matmul_gpu_dag
      module procedure matmul_gpu_dag_sp
      module procedure matmul_gpu_dag_dp
   end interface matmul_gpu_dag

   interface add_gpu
      module procedure add_gpu_sp
      module procedure add_gpu_dp
   end interface add_gpu

   interface add_cublas
      module procedure add_cublas_sp
      module procedure add_cublas_dp
   end interface add_cublas

   interface sum_gpu
      module procedure sum_gpu_sp
      module procedure sum_gpu_dp
   end interface sum_gpu

   interface inverse_gpu
      module procedure inverse_gpu_sp
      module procedure inverse_gpu_dp
   end interface inverse_gpu

   interface init_mat_gpu
      module procedure init_mat_gpu_sp
      module procedure init_mat_gpu_dp
   end interface init_mat_gpu

    interface zspectral_gpu
      module procedure zspectral_gpu_sp
      module procedure zspectral_gpu_dp
   end interface zspectral_gpu

 contains
   
 subroutine createGPU_sp(A)
   type(c_DNS ), intent(in) :: A     
   !$acc enter data create(A,A%val)
 end subroutine createGPU_sp
 
 subroutine createGPU_dp(A)
   type(z_DNS), intent(in) :: A     
   !$acc enter data create(A,A%val)
 end subroutine createGPU_dp 
 
 subroutine copyToGPU_sp(A)
   type(c_DNS ), intent(in) :: A     
   !$acc enter data copyin(A,A%val)
 end subroutine copyToGPU_sp

 subroutine copyToGPU_dp(A)
   type(z_DNS), intent(in) :: A     
   !$acc enter data copyin(A,A%val)
 end subroutine copyToGPU_dp
 
 subroutine copyFromGPU_sp(A)
   type(c_DNS ), intent(inout) :: A     
   !$acc update host(A%val)
 end subroutine copyFromGPU_sp

 subroutine copyFromGPU_dp(A)
   type(z_DNS), intent(inout) :: A     
   !$acc update host(A%val)
 end subroutine copyFromGPU_dp 
 
 subroutine deleteGPU_sp(A)
   type(c_DNS ), intent(inout) :: A     
   !$acc exit data delete(A%val,A)
 end subroutine deleteGPU_sp

 subroutine deleteGPU_dp(A)
   type(z_DNS), intent(inout) :: A     
   !$acc exit data delete(A%val,A)
 end subroutine deleteGPU_dp

 subroutine createAll_sp(Mat, nrow, ncol)
   type(c_DNS ) :: Mat 
   integer, intent(in) :: nrow, ncol

   call create(Mat,nrow,ncol)
   call createGPU(Mat)
   
 end subroutine createAll_sp

 subroutine createAll_dp(Mat, nrow, ncol)
   type(z_DNS) :: Mat 
   integer, intent(in) :: nrow, ncol

   call create(Mat,nrow,ncol)
   call createGPU(Mat)
   
 end subroutine createAll_dp

 subroutine destroyAll_sp(Mat)
   type(c_DNS ) :: Mat 
   call deleteGPU(Mat)
   call destroy(Mat)
 end subroutine destroyAll_sp

 subroutine destroyAll_dp(Mat)
   type(z_DNS) :: Mat 
   call deleteGPU(Mat)
   call destroy(Mat)
 end subroutine destroyAll_dp



 subroutine createMatGPU_sp(Mat)
   complex(sp), intent(in) :: Mat(:,:)      
   !$acc enter data create(Mat)
 end subroutine createMatGPU_sp

 subroutine createMatGPU_dp(Mat)
   complex(dp), intent(in) :: Mat(:,:)      
   !$acc enter data create(Mat)
 end subroutine createMatGPU_dp
 
 subroutine copyMatToGPU_sp(Mat)
  complex(sp), intent(in) :: Mat(:,:)      
  !$acc enter data copyin(Mat)
 end subroutine copyMatToGPU_sp

 subroutine copyMatToGPU_dp(Mat)
  complex(dp), intent(in) :: Mat(:,:)      
  !$acc enter data copyin(Mat)
 end subroutine copyMatToGPU_dp

 subroutine copyMatFromGPU_sp(Mat)
   complex(sp), intent(out) :: Mat(:,:)      
   !$acc update host(Mat)
 end subroutine copyMatFromGPU_sp

 subroutine copyMatFromGPU_dp(Mat)
   complex(dp), intent(out) :: Mat(:,:)      
   !$acc update host(Mat)
 end subroutine copyMatFromGPU_dp
 
 subroutine deleteMatGPU_sp(Mat)
   complex(sp), intent(in) :: Mat(:,:)      
   !$acc exit data delete(Mat)
 end subroutine deleteMatGPU_sp

 subroutine deleteMatGPU_dp(Mat)
   complex(dp), intent(in) :: Mat(:,:)      
   !$acc exit data delete(Mat)
 end subroutine deleteMatGPU_dp

 

 ! C = alpha*A*B + beta*C
 subroutine matmul_gpu_sp(hcublas,alpha,A,B,beta,C)
   type(cublasHandle), intent(in) :: hcublas
   complex(sp), intent(in) :: alpha      
   complex(sp), intent(in) :: A(:,:)      
   complex(sp), intent(in) :: B(:,:)      
   complex(sp), intent(in) :: beta     
   complex(sp), intent(inout) :: C(:,:)      

   integer :: istat, m,n,k
   m = size(A,1)
   k = size(B,1)
   n = size(C,2)
   !$acc data present(A, B, C)
   !$acc host_data use_device(A, B, C)
   istat = cublasCgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &
              & alpha, A, m, B, k, beta, C, m)
   !$acc end host_data
   !$acc end data

 end subroutine matmul_gpu_sp

 subroutine matmul_gpu_dp(hcublas,alpha,A,B,beta,C)
   type(cublasHandle), intent(in) :: hcublas
   complex(dp), intent(in) :: alpha      
   complex(dp), intent(in) :: A(:,:)      
   complex(dp), intent(in) :: B(:,:)      
   complex(dp), intent(in) :: beta     
   complex(dp), intent(inout) :: C(:,:)      

   integer :: istat, m,n,k
   m = size(A,1)
   k = size(B,1)
   n = size(C,2)
   !$acc data present(A, B, C)
   !$acc host_data use_device(A, B, C)
   istat = cublasZgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &
              & alpha, A, m, B, k, beta, C, m)
   !$acc end host_data
   !$acc end data

 end subroutine matmul_gpu_dp

 subroutine matmul_gpu_dag_sp(hcublas,alpha,A,B,beta,C)
 !subroutine matmul_sp_gpu_dag(hcublas,alpha,A,transA,B,transB,beta,C)
   type(cublasHandle), intent(in) :: hcublas
   complex(sp), intent(in) :: alpha      
   complex(sp), intent(in) :: A(:,:)      
   complex(sp), intent(in) :: B(:,:)
  ! type(cublasOperation), intent(in) :: transA, transB
   complex(sp), intent(in) :: beta     
   complex(sp), intent(inout) :: C(:,:)      

   integer :: istat, m,n,k
   m = size(A,1)
   k = size(B,1)
   n = size(C,2)
   !$acc data present(A, B, C)
   !$acc host_data use_device(A, B, C)
   istat = cublasCgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, &
              & alpha, A, m, B, k, beta, C, m)
   !$acc end host_data
   !$acc end data
 end subroutine matmul_gpu_dag_sp

 subroutine matmul_gpu_dag_dp(hcublas,alpha,A,B,beta,C)
   type(cublasHandle), intent(in) :: hcublas
   complex(dp), intent(in) :: alpha      
   complex(dp), intent(in) :: A(:,:)      
   complex(dp), intent(in) :: B(:,:)
  ! type(cublasOperation), intent(in) :: transA, transB
   complex(dp), intent(in) :: beta     
   complex(dp), intent(inout) :: C(:,:)      

   integer :: istat, m,n,k
   m = size(A,1)
   k = size(B,1)
   n = size(C,2)
   !$acc data present(A, B, C)
   !$acc host_data use_device(A, B, C)
   istat = cublasZgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, &
              & alpha, A, m, B, k, beta, C, m)
   !$acc end host_data
   !$acc end data
 end subroutine matmul_gpu_dag_dp

 
 !///////////////////////////////////////////////////////// 
 subroutine sum_gpu_sp(hcublas, Mat, summ)
   type(cublasHandle), intent(in) :: hcublas
   complex(sp), intent(in) :: Mat(:,:)
   real(sp), intent(out) :: summ
   !
   integer :: istat, n

   n=size(Mat,1)
 
   !$acc data present(Mat)
   !$acc host_data use_device(Mat)
   istat=cublasScasum(hcublas, n*n, Mat, 1, summ)
   !$acc end host_data
   !$acc end data
 end subroutine sum_gpu_sp

 subroutine sum_gpu_dp(hcublas, Mat, summ)
   type(cublasHandle), intent(in) :: hcublas
   complex(dp), intent(in) :: Mat(:,:)
   real(dp), intent(out) :: summ
   !
   integer :: istat, n

   n=size(Mat,1)
 
   !$acc data present(Mat)
   !$acc host_data use_device(Mat)
   istat=cublasDzasum(hcublas, n*n, Mat, 1, summ)
   !$acc end host_data
   !$acc end data
 end subroutine sum_gpu_dp

 subroutine copy_mat_gpu_sp(hcublas, A, B)
   type(cublasHandle), intent(in) :: hcublas
   complex(sp), intent(in) :: A(:,:)
   complex(sp), intent(out) :: B(:,:)

   integer :: istat, n 

   if (size(A) /= size(B)) then
      stop 'ERROR in copy_mat size mismatch'
   end if
   
   n=size(A)
   !$acc data present(A, B)
   !$acc host_data use_device(A,B)
   istat=cublasCcopy(hcublas, n, A, 1, B, 1)
   !$acc end host_data
   !$acc end data
 end subroutine copy_mat_gpu_sp

 subroutine copy_mat_gpu_dp(hcublas, A, B)
   type(cublasHandle), intent(in) :: hcublas
   complex(dp), intent(in) :: A(:,:)
   complex(dp), intent(out) :: B(:,:)

   integer :: istat, n 

   if (size(A) /= size(B)) then
      stop 'ERROR in copy_mat size mismatch'
   end if
   
   n=size(A)
   !$acc data present(A, B)
   !$acc host_data use_device(A,B)
   istat=cublasZcopy(hcublas, n, A, 1, B, 1)
   !$acc end host_data
   !$acc end data
 end subroutine copy_mat_gpu_dp

 subroutine inverse_gpu_sp(hcublas, hcusolver, A, Ainv, err)
   type(cublasHandle), intent(in) :: hcublas
   type(cusolverDnHandle), intent(in) :: hcusolver
   complex(sp), intent(in) :: A(:,:) 
   complex(sp), intent(inout) :: Ainv(:,:) 
   integer, intent(out) :: err

   complex(sp), allocatable :: LU(:,:) 
   complex(sp), allocatable, dimension(:) :: work 
   integer, allocatable, dimension(:) :: pivot
   integer :: n, ii, jj, lwork, istat, err1, err2
   !$acc declare present(A, Ainv)

   n= size(A,1)

   !$acc data present(A)
   istat=cusolverDnCgetrf_bufferSize(hcusolver,n,n,A,n,lwork)
   !$acc end data
   allocate(work(lwork))
   allocate(LU(n,n))
   allocate(pivot(n))
   
   !$acc kernels present(Ainv(n,n))
   Ainv = (0.0, 0.0)
   do jj = 1, n
     Ainv(jj,jj) = (1.0, 0.0)
   end do
   !$acc end kernels

   !$acc data create(LU, work, pivot, err1, err2) 
   !$acc host_data use_device(A, Ainv, LU, work, pivot, err1, err2)
   istat=cublasCcopy(hcublas, n*n, A, 1, LU, 1)
   istat=cusolverDnCgetrf(hcusolver,n,n,LU,n,work,pivot,err1)
   istat=cusolverDnCgetrs(hcusolver,CUBLAS_OP_N,n,n,LU,n,pivot,Ainv,n,err2)
   !$acc end host_data
   !$acc end data         
   err = err1

 end subroutine inverse_gpu_sp
 
 subroutine inverse_gpu_dp(hcublas, hcusolver, A, Ainv, err)
   type(cublasHandle), intent(in) :: hcublas
   type(cusolverDnHandle), intent(in) :: hcusolver
   complex(dp), intent(in) :: A(:,:) 
   complex(dp), intent(inout) :: Ainv(:,:) 
   integer, intent(out) :: err

   complex(dp), allocatable :: LU(:,:) 
   complex(dp), allocatable, dimension(:) :: work 
   integer, allocatable, dimension(:) :: pivot
   integer :: n, ii, jj, lwork, istat, err1, err2
   !$acc declare present(A, Ainv)

   n= size(A,1)

   !$acc data present(A)
   istat=cusolverDnZgetrf_bufferSize(hcusolver,n,n,A,n,lwork)
   !$acc end data
   allocate(work(lwork))
   allocate(LU(n,n))
   allocate(pivot(n))
   
   !$acc kernels present(Ainv(n,n))
   Ainv = (0.0, 0.0)
   do jj = 1, n
     Ainv(jj,jj) = (1.0, 0.0)
   end do
   !$acc end kernels

   !$acc data create(LU, work, pivot, err1, err2) 
   !$acc host_data use_device(A, Ainv, LU, work, pivot, err1, err2)
   istat=cublasZcopy(hcublas, n*n, A, 1, LU, 1)
   istat=cusolverDnZgetrf(hcusolver,n,n,LU,n,work,pivot,err1)
   istat=cusolverDnZgetrs(hcusolver,CUBLAS_OP_N,n,n,LU,n,pivot,Ainv,n,err2)
   !$acc end host_data
   !$acc end data         
   err = err1

 end subroutine inverse_gpu_dp

 ! -------------------------------------------------------------------------
 subroutine init_mat_gpu_sp(A)
   type(c_DNS ) :: A      

   integer :: jj

   !$acc kernels present(A)
   A%val = (0.0, 0.0)
   do jj = 1, size(A%val,1) 
     A%val(jj,jj) = (1.0, 0.0)
   end do
   !$acc end kernels

 end subroutine init_mat_gpu_sp

 subroutine init_mat_gpu_dp(A)
   type(z_DNS) :: A      

   integer :: jj

   !$acc kernels present(A)
   A%val = (0.0, 0.0)
   do jj = 1, size(A%val,1) 
     A%val(jj,jj) = (1.0, 0.0)
   end do
   !$acc end kernels

 end subroutine init_mat_gpu_dp

 ! -------------------------------------------------------------------------
 subroutine add_gpu_sp(C, A, B)
   type(c_DNS ), intent(in) :: A      
   type(c_DNS ), intent(in) :: B      
   type(c_DNS ), intent(inout) :: C      

   integer :: ii, jj, nrow, ncol
   nrow = size(A%val,1)
   ncol = size(A%val,2)

   !$acc kernels present(A, B, C)
   do jj = 1, ncol 
     do ii = 1, nrow 
        C%val(ii,jj) = A%val(ii,jj) + B%val(ii,jj)
     end do   
   end do
   !$acc end kernels
 
 end subroutine add_gpu_sp

 subroutine add_gpu_dp(C, A, B)
   type(z_DNS), intent(in) :: A      
   type(z_DNS), intent(in) :: B      
   type(z_DNS), intent(inout) :: C      

   integer :: ii, jj, nrow, ncol
   nrow = size(A%val,1)
   ncol = size(A%val,2)

   !$acc kernels present(A, B, C)
   do jj = 1, ncol 
     do ii = 1, nrow 
        C%val(ii,jj) = A%val(ii,jj) + B%val(ii,jj)
     end do   
   end do
   !$acc end kernels
 end subroutine add_gpu_dp

 subroutine zspectral_gpu_sp(Gr, A)
   type(c_DNS), intent(in) :: Gr            
   type(c_DNS), intent(inout) :: A      

   complex(sp) :: ii
   integer :: kk, jj, nrow, ncol
   nrow = size(A%val,1)
   ncol = size(A%val,2)

   ii = (0.0,1.0) 
   !$acc kernels present(Gr, A)
   do jj = 1, ncol 
      do kk = 1, nrow 
         A%val(kk,jj) = ii*(A%val(kk,jj) - conjg(A%val(jj,kk)))
      end do
   end do
   !$acc end kernels
 end subroutine zspectral_gpu_sp

 subroutine zspectral_gpu_dp(Gr, A)
   type(z_DNS), intent(in) :: Gr            
   type(z_DNS), intent(inout) :: A      

   complex(dp) :: ii
   integer :: kk, jj, nrow, ncol
   nrow = size(A%val,1)
   ncol = size(A%val,2)

   ii = (0.0,1.0) 
   !$acc kernels present(Gr, A)
   do jj = 1, ncol 
      do kk = 1, nrow 
         A%val(kk,jj) = ii*(A%val(kk,jj) - conjg(A%val(jj,kk)))
      end do
   end do
   !$acc end kernels
 end subroutine zspectral_gpu_dp

 subroutine add_cublas_sp(hcublas, C, A, alpha, B)
   !C = A + alpha B
   type(cublasHandle), intent(in) :: hcublas
   complex(sp), intent(in) :: A(:,:)      
   complex(sp), intent(in) :: B(:,:)      
   complex(sp), intent(inout) :: C(:,:)
   complex(sp), intent(in) :: alpha
   integer :: istat, n
   !complex(sp), parameter :: one = (1.0_sp,0.0_sp)
   if (size(A) /= size(B) .or. size(A) /= size(C)) then
      stop 'ERROR in add_cublas size mismatch'
   end if

   n = size(A)

   !$acc data present(A, B, C)
   !$acc host_data use_device(A, B, C)
   istat=cublasCcopy(hcublas, n, A, 1, C, 1)
   istat=cublasCaxpy(hcublas, n, alpha, B, 1, C, 1)
   !$acc end host_data
   !$acc end data

 end subroutine add_cublas_sp

 subroutine add_cublas_dp(hcublas, C, A, alpha, B)
   !C = A + alpha B
   type(cublasHandle), intent(in) :: hcublas
   complex(dp), intent(in) :: A(:,:)      
   complex(dp), intent(in) :: B(:,:)      
   complex(dp), intent(inout) :: C(:,:)
   complex(dp), intent(in) :: alpha
   integer :: istat, n
   !complex(sp), parameter :: one = (1.0_sp,0.0_sp)
   if (size(A) /= size(B) .or. size(A) /= size(C)) then
      stop 'ERROR in add_cublas size mismatch'
   end if

   n = size(A)

   !$acc data present(A, B, C)
   !$acc host_data use_device(A, B, C)
   istat=cublasZcopy(hcublas, n, A, 1, C, 1)
   istat=cublasZaxpy(hcublas, n, alpha, B, 1, C, 1)
   !$acc end host_data
   !$acc end data
 end subroutine add_cublas_dp

 subroutine copy_ESH_toGPU_sp(ESH)
   type(c_DNS), dimension(:,:), intent(in) :: ESH
   integer :: ii, nbl

   nbl = size(ESH,1)
   call copyToGPU(ESH(1,1))
   do ii=1,nbl
      call copyToGPU(ESH(ii,ii))
      call copyToGPU(ESH(ii-1,ii))
      call copyToGPU(ESH(ii,ii-1))
   end do
end subroutine copy_ESH_toGPU_sp 

subroutine copy_ESH_toGPU_dp(ESH)
   type(z_DNS), dimension(:,:), intent(in) :: ESH
   integer :: ii, nbl

   nbl = size(ESH,1)
   call copyToGPU(ESH(1,1))
   do ii=1,nbl
      call copyToGPU(ESH(ii,ii))
      call copyToGPU(ESH(ii-1,ii))
      call copyToGPU(ESH(ii,ii-1))
   end do
end subroutine copy_ESH_toGPU_dp

end module cudautils
