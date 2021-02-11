#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#ifdef __INTEL_COMPILER
  #define MKL_Complex16 double complex
  #include "mkl.h"
#endif
#include <time.h>
#include <omp.h>

//#include "global_parameters.h"
//#include "simulation_parameters.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define KK(im, iK, iQ) ( KK[ im + Np * iK + Np * NK * iQ ] )


// =========================================================
// self-energy for polar optical electron-phonon interaction
// =========================================================

int self_energy(
  // communicator of the 2d cartesian grid  
  MPI_Comm comm2d,
  // number of rows of the matrix block (Np*Mp) 
  int Np, 
  // number of cols of the matrix block (Np*Mp) 
  int Mp, 
  // global number of K-points
  int NK,
  // global number of E-points
  int NE,
  // local number of K-points
  int NKloc,
  // local number of E-points
  int NEloc,
  // hbar wq in energy grid units 
  int iEhbaromegaLO,
  // array of pointers to G (double complex *)
  void **GG,
  // array of pointers to Sigma (double complex *)
  void **Sigma,
  // a series of buffers size Np*Np (can they be reused ?) 
  double complex *sbuff1,
  double complex *sbuff2,
  double complex *rbuff1, 
  double complex *rbuff2,
  double complex *sbuffH, 
  double complex *rbuffH,
  // fac1 * G(E-wq)
  double fac1,
  // fac2 * G(E+wq)
  double fac2,
  // Pretubulated array KK(iQ, iK, |zi-zj|) 
  double *KK,
)
{
  int myid, size;
  int iK,iE, iQ;
  int iKglo,iEglo,iQglo, iKglo2;
  int iz,jz, im, in;
  int coords[2], coordsH[2], dims[2];
  int iEminus, iEplus;
  
  int msource, mdest;
  int psource, pdest;
  int hsource, hdest;
  int ndiff;
  int ndims;

  double complex *pbuff1, *pbuff2, *pGG, *pSigma;

  // dimensions of the cartesian grid
  ndims = 2; 
  dims[0] = NK/NKloc;
  dims[1] = NE/NEloc;

  MPI_Request rqE[4];
  MPI_Status statusE[4];
  
  MPI_Request rqH[2];
  MPI_Status statusH[2];
  
  MPI_Comm_size(comm2d, &size );
  MPI_Comm_rank(comm2d, &myid);
  MPI_Cart_coords(comm2d,myid,ndims,coords);
  
  // CAREFUL! The meaning of K and Q is reversed compared to the notes!
  //
  // Sigma_ij(iQ, iE) = Sum_iK   F(iQ, iK, |z_i-z_j|) * 
  //                        * (fac1 * GG_ij(iK, E-wq) + fac2 * GG_ij(iK, E+wq))
  //
  for( iK=0; iK<NKloc; iK++ )
  {
    iKglo=iK+coords[0]*NKloc;
    for( iE=0; iE<NEloc; iE++ )
    {
      iEglo=iE+coords[1]*NEloc;
      
      ////////////////////////////////////////////////////////////////////////////////////
      // Communications of G(k, E-hwq) 
      ////////////////////////////////////////////////////////////////////////////////////

      // checks if iE-iEhbaromegaLO is on the same processor => no communication
      if( (((int)iE - iEhbaromegaLO) >=0 && (iE - iEhbaromegaLO) < NEloc) || iEglo < NEloc)
      {
        iEminus=iE-iEhbaromegaLO;
        if(iEminus<=0) {iEminus = 0;}
        // pbuff1 points to G(k,E-wq)
        //pbuff1 = &Gless[(iK*NEloc+iEminus)*Np*Np]  -> GG[iK,iE]
        pbuff1 = (double complex *) GG[iK*NEloc+iEminus];
        mdest=MPI_PROC_NULL;
        msource=MPI_PROC_NULL;
      }
      
      if(dims[1] > 1)
      {
        // get global iEminus
        iEminus = iEglo-iEhbaromegaLO;
        if( iEhbaromegaLO >= NEloc )
          ndiff = (iEglo-iEminus)/ NEloc;
        else
          ndiff = (NEloc+iEglo-iEminus)/ NEloc;
      
        // get destination process 
        MPI_Cart_shift( comm2d, 1, ndiff, &msource, &mdest );
        
        if(mdest != MPI_PROC_NULL && iE<iEhbaromegaLO) 
        {
          // get local iEmins in destination process and send there
          iEminus = (NEloc+iE-iEhbaromegaLO) % NEloc;
          pGG = (double complex *) GG[iK*NEloc+iEminus];
          MPI_Isend(pGG, Mp*Np, MPI_DOUBLE_COMPLEX, mdest, 41, comm2d, &rqE[0]);
        }

        // corresponding recv        
        if(msource != MPI_PROC_NULL && iE<iEhbaromegaLO) 
        {
          MPI_Irecv(rbuff1, Mp*Np, MPI_DOUBLE_COMPLEX, msource, 41, comm2d, &rqE[1]);
          // pbuff1 points to received buffer with G(k,E-wq)
          pbuff1 = rbuff1;
        }
      }
      
      ////////////////////////////////////////////////////////////////////////////////////
      // Communications of G(k, E+hwq) 
      ////////////////////////////////////////////////////////////////////////////////////
      
      // checks if iE+iEhbaromegaLO is on the same processor => no communication
      if( (((int)iE + iEhbaromegaLO) >=0 && (iE + iEhbaromegaLO) < NEloc) || iEglo >= NE-NEloc)
      {
        iEplus = iE+iEhbaromegaLO;
        if(iEglo+iEhbaromegaLO>=NE) {iEplus = NEloc-1;}
        //pbuff2 = &Gless[(iK*NEloc+iEplus)*Np*Np];
        pbuff2 = (double complex *) GG[iK*NEloc+iEplus];
        pdest=MPI_PROC_NULL;
        psource=MPI_PROC_NULL;
      }
      
      if(dims[1] > 1)
      {
        iEplus = iEglo + iEhbaromegaLO;
        if( iEhbaromegaLO >= NEloc )
          ndiff = (iEplus-iEglo) / NEloc;
        else
          ndiff = (iEplus+NEloc-iEglo)/ NEloc;

        MPI_Cart_shift( comm2d, 1, -ndiff, &psource, &pdest );
        // get local iEmins in destination process and send there
        if(pdest != MPI_PROC_NULL && iE >= NEloc-iEhbaromegaLO)
        {
          iEplus = iE - (NEloc - iEhbaromegaLO);
          pGG = (double complex *) GG[iK*NEloc+iEplus];
          MPI_Isend(pGG, Mp*Np, MPI_DOUBLE_COMPLEX, pdest, 42, comm2d, &rqE[2]);
        }
        // recv from source
        if(psource != MPI_PROC_NULL && iE >= NEloc-iEhbaromegaLO)
        {
          MPI_Irecv(rbuff2, Mp*Np,MPI_DOUBLE_COMPLEX,psource,42,comm2d,&rqE[3]);
          // pbuff2 points G(k,E+hwq)
          pbuff2 = rbuff2;
        }
        
      }
      
      if(dims[1] > 1)
      {
        if(mdest != MPI_PROC_NULL  && iE<iEhbaromegaLO) MPI_Wait(&rqE[0],&statusE[0]);
        if(msource != MPI_PROC_NULL  && iE<iEhbaromegaLO) MPI_Wait(&rqE[1],&statusE[1]);
        
        if(pdest != MPI_PROC_NULL && iE >= NEloc-iEhbaromegaLO) MPI_Wait(&rqE[2],&statusE[2]);
        if(psource != MPI_PROC_NULL && iE >= NEloc-iEhbaromegaLO) MPI_Wait(&rqE[3],&statusE[3]);
      }
      
      ////////////////////////////////////////////////////////////////////////////////////
      // Communications of k-points 
      ////////////////////////////////////////////////////////////////////////////////////
      
      // Update local
      #pragma omp parallel for private(iQ,jz,iz,iQglo,im) collapse(2)
      for( iQ=0; iQ<NKloc; iQ++ )
      {
        iQglo=iQ+coords[0]*NKloc;
        for( jz=0; jz<Mp; jz++ )
        {
          for( iz=0; iz<Np; iz++ )
          {
            im = abs(iz-jz);
            pSigma = (double complex *) Sigma[iQ*NEloc + iE];
            pSigma[jz*Np + iz] += KK(im, iKglo, iQ) * 
                                 (fac1 * pbuff1[jz*Np+iz] + fac2 * pbuff2[jz*Np+iz]); 
            
            // F(iQ, iK, |iz-jz|) = F[Np*NK*iQglo + Np*iKglo + im) 
            //sigless[(iQ*NEloc+iE)*Np*Np+jz*Np+iz] += dK * qe * hbaromegaLO/(4.0*pow(M_PI,2.0)) * 
            //  (1.0/(epsinfr*eps0) - 1.0/(eps0r*eps0)) * K[iKglo] * F[ im+ Np*iKglo + Np*NK*iQglo ] * 
            //  Mtilde * (fac1 * pbuff1[jz*Np+iz] + fac2 * pbuff2[jz*Np+iz]); 
          }
        }
      }
      
      
      if(dims[0] > 1)
      {
        //printf("Communication for Q integration\n");
        #pragma omp parallel for private(iz,jz)
        for( jz=0; jz<Mp; jz++ )
          for( iz=0; iz<Np; iz++ )
            sbuffH[jz*Np+iz] = fac1 * pbuff1[jz*Np+iz] + fac2 * pbuff2[jz*Np+iz];
      
      
        for( in=1; in<dims[0]; in++ )
        {
          MPI_Cart_shift(comm2d, 0, in, &hsource, &hdest);
          
          MPI_Cart_coords(comm2d, hsource, ndims, coordsH);
         
          iKglo2=iK+coordsH[0]*NKloc;
          
          // send to in+1%dims[0] = hdest
          MPI_Isend(sbuffH, Mp*Np, MPI_DOUBLE_COMPLEX, hdest, 43, comm2d, &rqH[0]);
          // recv from in+dims[0]%dims[0] = hsource
          MPI_Irecv(rbuffH, Mp*Np, MPI_DOUBLE_COMPLEX, hsource, 43, comm2d, &rqH[1]);
          
          MPI_Wait(&rqH[1], &statusH[1]);
          
          #pragma omp parallel for private(iQ,jz,iz,iQglo,im) collapse(2)
          for( iQ=0; iQ<NKloc; iQ++ )
          {
            iQglo=iQ+coords[0]*NKloc;
            for( jz=0; jz<Mp; jz++ )
            {
              for( iz=0; iz<Np; iz++ )
              {
                im = abs(iz-jz);
                pSigma = (double complex *) Sigma[iQ*NEloc + iE];
                pSigma[jz*Np + iz] += KK(im, iKglo2, iQ) * rbuffH[jz*Np+iz]; 
                //sigless[(iQ*NEloc+iE)*Np*Np+jz*Np+iz] += dK * qe * hbaromegaLO/(4.0*pow(M_PI,2.0)) * 
                //  (1.0/(epsinfr*eps0) - 1.0/(eps0r*eps0)) * K[iKglo2] * F[ im+ Np*iKglo2 + Np*NK*iQglo ] * 
                //  Mtilde * rbuffH[jz*Np+iz];
              }
            }
          }
          
          MPI_Wait(&rqH[0],&statusH[0]);
        }
      }
      
    }
  }
  
  MPI_Barrier(comm2d);
  
  return 0;

}
