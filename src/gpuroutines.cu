#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "cusolverDn.h"

#define BLOCK_SIZE 1024
#define TILE_DIM 32 
#define BLOCK_ROWS 8

__global__ void CaddKernel(cuComplex* c, const cuComplex alpha, const cuComplex* a, const cuComplex beta, const cuComplex* b, int size)
{
   int i = blockIdx.x * blockDim.x + threadIdx.x;
   if (i < size) {
      c[i].x = (alpha.x * a[i].x - alpha.y * a[i].y) + (beta.x * b[i].x - beta.y * b[i].y);
      c[i].y = (alpha.x * a[i].y + alpha.y * a[i].x) + (beta.x * b[i].y + beta.y * b[i].x);
   }
}

__global__ void ZaddKernel(cuDoubleComplex* c, const cuDoubleComplex alpha, const cuDoubleComplex* a, const cuDoubleComplex beta, const cuDoubleComplex* b, int size)
{
   int i = blockIdx.x * blockDim.x + threadIdx.x;
   if (i < size) {
      c[i].x = (alpha.x * a[i].x - alpha.y * a[i].y) + (beta.x * b[i].x - beta.y * b[i].y);
      c[i].y = (alpha.x * a[i].y + alpha.y * a[i].x) + (beta.x * b[i].y + beta.y * b[i].x);
   }
}

/*
__global__ void hermitian(cuComplex *odata, const cuComplex *idata)
{
  __shared__ cuComplex tile[TILE_DIM][TILE_DIM];
  
  int x = blockIdx.x * TILE_DIM + threadIdx.x;
  int y = blockIdx.y * TILE_DIM + threadIdx.y;
  int width = gridDim.x * TILE_DIM;

  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
     tile[threadIdx.y+j][threadIdx.x] = idata[(y+j)*width + x];

  __syncthreads();

  x = blockIdx.y * TILE_DIM + threadIdx.x;  // transpose block offset
  y = blockIdx.x * TILE_DIM + threadIdx.y;

  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
  {	   
     odata[(y+j)*width + x].x = tile[threadIdx.x][threadIdx.y + j].x;
     odata[(y+j)*width + x].y = -tile[threadIdx.x][threadIdx.y + j].y;
  }
}
*/
__global__ void CinitKernel(cuComplex *a, int nrow) {
    
    int size;	 
    int i = blockDim.x*blockIdx.x + threadIdx.x;

    size = nrow*nrow;
    if(i < size) {
          if(i%(nrow+1) == 0){
              a[i].x = 1.0;
              a[i].y = 0.0;
	     }
          else{
              a[i].x = 0.0;
              a[i].y = 0.0;
	   }
    }
}

__global__ void ZinitKernel(cuDoubleComplex *a, int nrow) {
    
    int size;	 
    int i = blockDim.x*blockIdx.x + threadIdx.x;

    size = nrow*nrow;
    if(i < size) {
          if(i%(nrow+1) == 0){
              a[i].x = 1.0;
              a[i].y = 0.0;
	     }
          else{
              a[i].x = 0.0;
              a[i].y = 0.0;
	   }
    }
}

__global__ void CtraceKernel(cuComplex *a, int nrow, cuComplex *trace) {
    
    int size;	 
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    
    size = nrow*nrow;
    if(i < size) {
          if(i%(nrow+1) == 0){
             trace[i].x = a[i].x;
             trace[i].y = 0.0;
             }
          else{
             trace[i].x = 0.0;
             trace[i].y = 0.0;
	     }
    }
}

__global__ void ZtraceKernel(cuDoubleComplex *a, int nrow, cuDoubleComplex *trace) {
    
    int size;	 
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    
    size = nrow*nrow;
    if(i < size) {
          if(i%(nrow+1) == 0){
             trace[i].x = a[i].x;
             trace[i].y = 0.0;
             }
          else{
             trace[i].x = 0.0;
             trace[i].y = 0.0;
	     }
    }
}
/*~-~-~-~-~-~-~-~-~-~-~-~-~-~ DATA MOVEMENT  ROUTINES -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-*/

extern "C" int cu_createMat(void **d_A, int bytecount)
{
  cudaError_t err;
  err = cudaMalloc(d_A, bytecount);
  printf("GPU Address: %p \n",*d_A);
  return err;
}

extern "C" int cu_copyMatH2D(void *h_A, void *d_A, int bytecount)
{
  cudaError_t err;
  printf("copy %p to %p\n",h_A,d_A);
  err = cudaMemcpy(d_A, h_A, bytecount, cudaMemcpyHostToDevice);  
  return err;
}	

extern "C" int cu_copyMatD2H(void *h_A, void *d_A, int bytecount)
{
  cudaError_t err;
  printf("copy %p to %p\n",d_A,h_A);
  err = cudaMemcpy(h_A, d_A, bytecount, cudaMemcpyDeviceToHost);  
  return err;
}	

extern "C" int cu_deleteMat(void *d_A)
{
  cudaError_t err;
  printf("add_free: %p",d_A);
  err = cudaFree(d_A);  
  return err;
}	

/*~-~-~-~-~-~-~-~-~-~-~-~-~-~ INIT/FINAL ROUTINES -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-*/

extern "C" int cu_cublasInit(cublasHandle_t *hcublas)
{
  cublasStatus_t err;
  err = cublasCreate(hcublas);
  if (err != 0){
    printf("cublas create error: %d\n",err);
  }
  printf("hcublas Addr: %p \n",*hcublas);
  return err;
}

extern "C" int cu_cublasFinalize(cublasHandle_t hcublas)
{
  cublasStatus_t err;
  err = cublasDestroy(hcublas); 
  return err;
}

extern "C" int cu_cusolverInit(cusolverDnHandle_t *hcusolver)
{
  cusolverStatus_t err;
  err = cusolverDnCreate(hcusolver);
  if (err != 0){
    printf("cusolver create error: %d\n",err);
  }
  printf("hcusolver Addr: %p \n",*hcusolver);
  return err;
}

extern "C" int cu_cusolverFinalize(cusolverDnHandle_t hcusolver)
{
  cusolverStatus_t err;
  err = cusolverDnDestroy(hcusolver); 
  return err;
}

/*~-~-~-~-~-~-~-~-~-~-~-~-~-~ MATRIX ROUTINES -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-*/

extern "C" int cu_CmultMat(cublasHandle_t hcublas, int m, int n, int k, cuComplex *alpha, void *d_A, void *d_B, cuComplex *beta, void *d_C, int dagger)
{
  cuComplex *pdA, *pdB, *pdC; 	

  printf("A: %p B: %p C: %p\n",d_A,d_B,d_C);
  pdA=(cuComplex *) d_A;
  pdB=(cuComplex *) d_B;
  pdC=(cuComplex *) d_C;
  cublasStatus_t err;
  if (dagger == 0){
     err = cublasCgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, pdA, m, pdB, k, beta, pdC, m);
  }
  if (dagger == 1){
     err = cublasCgemm(hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, k, alpha, pdA, k, pdB, k, beta, pdC, m);
  }
  if (dagger == 2){  
     err = cublasCgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, alpha, pdA, m, pdB, n, beta, pdC, m);
  }
  return err;
}

extern "C" int cu_ZmultMat(cublasHandle_t hcublas, int m, int n, int k, cuDoubleComplex *alpha, void *d_A, void *d_B, cuDoubleComplex *beta, void *d_C, int dagger)
{
  cuDoubleComplex *pdA, *pdB, *pdC; 	

  printf("A: %p B: %p C: %p\n",d_A,d_B,d_C);
  pdA=(cuDoubleComplex *) d_A;
  pdB=(cuDoubleComplex *) d_B;
  pdC=(cuDoubleComplex *) d_C;
  cublasStatus_t err;
  if (dagger == 0){
     err = cublasZgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, pdA, m, pdB, k, beta, pdC, m);
  }
  if (dagger == 1){
     err = cublasZgemm(hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, k, alpha, pdA, k, pdB, k, beta, pdC, m);
  }
  if (dagger == 2){  
     err = cublasZgemm(hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, alpha, pdA, m, pdB, n, beta, pdC, m);
  }
  return err;
}

extern "C" int cu_Cinverse(cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void *d_A, void *d_Ainv, int N)
{
   cudaError_t cudaStatus;
   cusolverStatus_t cusolverStatus;
   cublasStatus_t cublasStatus;
   // declare arrays on the device
   cuComplex  *pdA , *pdAinv, *d_LU, *d_Work; 
  
   pdA = (cuComplex *) d_A;
   pdAinv = (cuComplex *) d_Ainv;
   // coeff . matrix , rhs , workspace
   int *d_pivot , *d_info , Lwork ; // pivots , info , worksp . size
   int info_gpu = 0;

   // compute buffer size and prep . memory
   cusolverStatus = cusolverDnCgetrf_bufferSize( hcusolver, N , N , pdA , N , &Lwork);
   // prepare memory on the device

   cudaStatus = cudaMalloc(( void **)& d_LU, N*N*sizeof(cuComplex));
   cudaStatus = cudaMalloc(( void **)& d_pivot , N*sizeof(int));
   cudaStatus = cudaMalloc(( void **)& d_info , sizeof(int));
   // copy d_LU <- pdA
   cublasStatus = cublasCcopy(hcublas, N*N, pdA, 1, d_LU, 1);

   cudaStatus = cudaMalloc(( void **)& d_Work , Lwork*sizeof(cuComplex));
  
   // LU factorization of d_A , with partial pivoting and row
   // interchanges ; row i is interchanged with row d_pivot ( i );
   cusolverStatus = cusolverDnCgetrf(hcusolver, N, N, d_LU, N, d_Work, d_pivot, d_info);
  
   // use the LU factorization to solve the system d_LU * x = d_Ainv ;
   // the solution overwrites d_Ainv
   cusolverStatus = cusolverDnCgetrs(hcusolver, CUBLAS_OP_N, N, N, d_LU, N, d_pivot, pdAinv, N, d_info);

   cudaStatus = cudaMemcpy(&info_gpu , d_info , sizeof(int), cudaMemcpyDeviceToHost);
   // d_info -> info_gpu
   cudaStatus = cudaFree(d_pivot);
   cudaStatus = cudaFree(d_info);
   cudaStatus = cudaFree(d_Work);
   cudaStatus = cudaFree(d_LU);
   return cudaStatus;
}

extern "C" int cu_Zinverse(cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void *d_A, void *d_Ainv, int N)
{
   cudaError_t cudaStatus;
   cusolverStatus_t cusolverStatus;
   cublasStatus_t cublasStatus; 
   // declare arrays on the device
   cuDoubleComplex  *pdA , *pdAinv, *d_LU, *d_Work; 
  
   pdA = (cuDoubleComplex *) d_A;
   pdAinv = (cuDoubleComplex *) d_Ainv;
   // coeff . matrix , rhs , workspace
   int *d_pivot , *d_info , Lwork ; // pivots , info , worksp . size
   int info_gpu = 0;

   // compute buffer size and prep . memory
   cusolverStatus = cusolverDnZgetrf_bufferSize( hcusolver, N , N , pdA , N , &Lwork);
   // prepare memory on the device

   cudaStatus = cudaMalloc(( void **)& d_LU, N*N*sizeof(cuDoubleComplex));
   cudaStatus = cudaMalloc(( void **)& d_pivot , N*sizeof(int));
   cudaStatus = cudaMalloc(( void **)& d_info , sizeof(int));
   // copy d_LU <- pdA
   cublasStatus = cublasZcopy(hcublas, N*N, pdA, 1, d_LU, 1);

   cudaStatus = cudaMalloc(( void **)& d_Work , Lwork*sizeof(cuDoubleComplex));
  
   // LU factorization of d_A , with partial pivoting and row
   // interchanges ; row i is interchanged with row d_pivot ( i );
   cusolverStatus = cusolverDnZgetrf(hcusolver, N, N, d_LU, N, d_Work, d_pivot, d_info);
  
   // use the LU factorization to solve the system d_LU * x = d_Ainv ;
   // the solution overwrites d_Ainv
   cusolverStatus = cusolverDnZgetrs(hcusolver, CUBLAS_OP_N, N, N, d_LU, N, d_pivot, pdAinv, N, d_info);

   cudaStatus = cudaMemcpy(&info_gpu , d_info , sizeof(int), cudaMemcpyDeviceToHost);
   // d_info -> info_gpu
   cudaStatus = cudaFree(d_pivot);
   cudaStatus = cudaFree(d_info);
   cudaStatus = cudaFree(d_Work);
   cudaStatus = cudaFree(d_LU);
   return cudaStatus;
}

extern "C" int cu_Ckernelsum(void *d_C, cuComplex *alpha, void *d_A, cuComplex *beta, void *d_B, int size)
{
   int NumBlocks;
   cuComplex *pdA = (cuComplex *) d_A;
   cuComplex *pdB = (cuComplex *) d_B;
   cuComplex *pdC = (cuComplex *) d_C;

   NumBlocks = (size/BLOCK_SIZE)+1;

   CaddKernel<<<NumBlocks,BLOCK_SIZE>>>(pdC, *alpha, pdA, *beta, pdB, size);

   return 0; 
}

extern "C" int cu_Zkernelsum(void *d_C, cuDoubleComplex *alpha, void *d_A, cuDoubleComplex *beta, void *d_B, int size)
{
   int NumBlocks;
   cuDoubleComplex *pdA = (cuDoubleComplex *) d_A;
   cuDoubleComplex *pdB = (cuDoubleComplex *) d_B;
   cuDoubleComplex *pdC = (cuDoubleComplex *) d_C;

   NumBlocks = (size/BLOCK_SIZE)+1;

   ZaddKernel<<<NumBlocks,BLOCK_SIZE>>>(pdC, *alpha, pdA, *beta, pdB, size);

   return 0; 
}

extern "C" int cu_Cmatsum(cublasHandle_t hcublas, int m, int n, cuComplex *alpha, void *d_A, cuComplex *beta, void *d_B, void *d_C, int dagger)
{
   //m number of rows of matrix op(A) and C
   //n number of columns of matrix op(B) and C  	
   cuComplex *pdA = (cuComplex *) d_A;
   cuComplex *pdB = (cuComplex *) d_B;
   cuComplex *pdC = (cuComplex *) d_C;

   cublasStatus_t err;
   if (dagger == 0) {
      err = cublasCgeam(hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, alpha, pdA, m, beta, pdB, m, pdC, m);
      }
   if (dagger == 1) {
      err = cublasCgeam(hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, alpha, pdA, n, beta, pdB, m, pdC, m);
      }
   if (dagger == 2) {
      err = cublasCgeam(hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, alpha, pdA, m, beta, pdB, n, pdC, m);
      }
   return err; 
}

extern "C" int cu_Zmatsum(cublasHandle_t hcublas, int m, int n, cuDoubleComplex *alpha, void *d_A, cuDoubleComplex *beta, void *d_B, void *d_C, int dagger)
{
   //m number of rows of matrix op(A) and C
   //n number of columns of matrix op(B) and C  	
   cuDoubleComplex *pdA = (cuDoubleComplex *) d_A;
   cuDoubleComplex *pdB = (cuDoubleComplex *) d_B;
   cuDoubleComplex *pdC = (cuDoubleComplex *) d_C;

   cublasStatus_t err;
   if (dagger == 0) {
      err = cublasZgeam(hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, alpha, pdA, m, beta, pdB, m, pdC, m);
      }
   if (dagger == 1) {
      err = cublasZgeam(hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, alpha, pdA, n, beta, pdB, m, pdC, m);
      }
   if (dagger == 2) {
      err = cublasZgeam(hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, alpha, pdA, m, beta, pdB, n, pdC, m);
      }
   return err; 
}

extern "C" int cu_Cinitmat( void *d_A, int nrow)
{
   int NumBlocks;
   int size = nrow*nrow;
   cuComplex *pdA = (cuComplex *) d_A;

   NumBlocks = (size/BLOCK_SIZE)+1;

   CinitKernel<<<NumBlocks,BLOCK_SIZE>>>(pdA, nrow);

   return 0; 
}

extern "C" int cu_Zinitmat( void *d_A, int nrow)
{
   int NumBlocks;
   int size = nrow*nrow;
   cuDoubleComplex *pdA = (cuDoubleComplex *) d_A;

   NumBlocks = (size/BLOCK_SIZE)+1;

   ZinitKernel<<<NumBlocks,BLOCK_SIZE>>>(pdA, nrow);

   return 0; 
}

extern "C" float cu_Ctrace(cublasHandle_t hcublas, void *d_A, int nrow)
{
   cudaError_t cudaStatus;
   cublasStatus_t err;
   
   int NumBlocks;
   float result;
   int size = nrow*nrow;
   cuComplex *d_work;
   cuComplex *pdA = (cuComplex *) d_A;

   NumBlocks = (size/BLOCK_SIZE)+1;
   
   cudaStatus = cudaMalloc(( void **)& d_work, size*sizeof(cuComplex));

   CtraceKernel<<<NumBlocks,BLOCK_SIZE>>>(pdA, nrow, d_work);
   err = cublasScasum(hcublas, size, d_work, 1, &result);

   cudaStatus = cudaFree(d_work);
   
   return result;
}

extern "C" double cu_Ztrace(cublasHandle_t hcublas, void *d_A, int nrow)
{
   cudaError_t cudaStatus;
   cublasStatus_t err;
   
   int NumBlocks;
   double result;
   int size = nrow*nrow;
   cuDoubleComplex *d_work;
   cuDoubleComplex *pdA = (cuDoubleComplex *) d_A;

   NumBlocks = (size/BLOCK_SIZE)+1;
   
   cudaStatus = cudaMalloc(( void **)& d_work, size*sizeof(cuDoubleComplex));

   ZtraceKernel<<<NumBlocks,BLOCK_SIZE>>>(pdA, nrow, d_work);
   err = cublasDzasum(hcublas, size, d_work, 1, &result);

   cudaStatus = cudaFree(d_work);
   
   return result;
}

extern "C" int cu_Cmatcopy(cublasHandle_t hcublas,  void *d_A,  void *d_B, int N)
{
   //m number of rows of matrix op(A) and C
   //n number of columns of matrix op(B) and C  	
   cuComplex *pdA = (cuComplex *) d_A;
   cuComplex *pdB = (cuComplex *) d_B;

   cublasStatus_t err;
   
   err = cublasCcopy(hcublas, N*N, pdA, 1, pdB, 1);
   return err; 
}

extern "C" int cu_Zmatcopy(cublasHandle_t hcublas,  void *d_A,  void *d_B, int size)
{
   //m number of rows of matrix op(A) and C
   //n number of columns of matrix op(B) and C  	
   cuDoubleComplex *pdA = (cuDoubleComplex *) d_A;
   cuDoubleComplex *pdB = (cuDoubleComplex *) d_B;

   cublasStatus_t err;
   
   err = cublasZcopy(hcublas, size, pdA, 1, pdB, 1);
   return err; 
}

