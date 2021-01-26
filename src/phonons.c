#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#ifdef __INTEL_COMPILER
  #define MKL_Complex16 double complex
  #include "mkl.h"
#endif
#include <time.h>
#include <omp.h>

#include "global_parameters.h"
#include "simulation_parameters.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// =========================================================
// self-energy for polar optical electron-phonon interaction
// =========================================================

int self_energy_polar_optical_phonon_scattering(
  int Np,
  int NK,
  int NE,
  int NKloc,
  int NEloc,
  int iEhbaromegaLO,
  double complex *GG,
  double complex *sigma,
  double complex *sbuff1,
  double complex *sbuff2,
  double complex *rbuff1, 
  double complex *rbuff2,
  double complex *sbuffH, 
  double complex *rbuffH,
  MPI_Comm comm2d, 
  MPI_Comm comm2d_per,
  double fac1,
  double fac2
)
{
  int myid, size;
  int iK,iE, iQ;
  int iKglo,iEglo,iQglo, iKglo2;
  int iz,jz, im, in;
  int coords[2], coordsH[2];
  int iEminus, iEplus;
  
  int msource, mdest;
  int psource, pdest;
  int hsource, hdest;
  int ndiff;
  
  MPI_Request rqE[4];
  MPI_Status statusE[4];
  
  MPI_Request rqH[2];
  MPI_Status statusH[2];
  
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &myid);
  MPI_Cart_coords(comm2d,myid,ndims,coords);
  
  
  for( iK=0; iK<NKloc; iK++ )
  {
    iKglo=iK+coords[0]*NKloc;
    for( iE=0; iE<NEloc; iE++ )
    {
      iEglo=iE+coords[1]*NEloc;
      
      
      if( (((int)iE - iEhbaromegaLO) >=0 && (iE - iEhbaromegaLO) < NEloc) || iEglo < NEloc)
      {
        // No communication neccassary for iEminus
        iEminus=iE-iEhbaromegaLO;
        if(iEminus<=0) {iEminus = 0;}
        pbuff1 = &Gless[(iK*NEloc+iEminus)*Np*Np];
        mdest=MPI_PROC_NULL;
        msource=MPI_PROC_NULL;
      }
      
      if(dims[1] > 1)
      {
        iEminus = iEglo-iEhbaromegaLO;
        if( iEhbaromegaLO >= NEloc )
          ndiff = (iEglo-iEminus)/ NEloc;
        else
          ndiff = (NEloc+iEglo-iEminus)/ NEloc;
        MPI_Cart_shift( comm2d, 1, ndiff, &msource, &mdest );
        // send to dest
        if(mdest != MPI_PROC_NULL && iE<iEhbaromegaLO) 
        {
          iEminus = (NEloc+iE-iEhbaromegaLO) % NEloc;
          MPI_Isend(&Gless[(iK*NEloc+iEminus)*Np*Np],Np*Np,MPI_DOUBLE_COMPLEX,mdest,41,comm2d,&rqE[0]);
        }
        // recv from source
        if(msource != MPI_PROC_NULL && iE<iEhbaromegaLO) 
        {
          MPI_Irecv(rbuff1,Np*Np,MPI_DOUBLE_COMPLEX,msource,41,comm2d,&rqE[1]);
          // pointer
          pbuff1 = rbuff1;
        }
      }
      
      
      
      
      
      
      
      
      
      if( (((int)iE + iEhbaromegaLO) >=0 && (iE + iEhbaromegaLO) < NEloc) || iEglo >= NE-NEloc)
      {
        // No communication neccassary for iEplus
        iEplus = iE+iEhbaromegaLO;
        if(iEglo+iEhbaromegaLO>=NE) {iEplus = NEloc-1;}
        pbuff2 = &Gless[(iK*NEloc+iEplus)*Np*Np];
        pdest=MPI_PROC_NULL;
        psource=MPI_PROC_NULL;
      }
      
      if(dims[1] > 1)
      {
        iEplus = iEglo + iEhbaromegaLO;
        if( iEhbaromegaLO >= NEloc )
        {
          ndiff = (iEplus-iEglo) / NEloc;
        }
        else
        {
          ndiff = (iEplus+NEloc-iEglo)/ NEloc;
        }
        MPI_Cart_shift( comm2d, 1, -ndiff, &psource, &pdest );
        // send to dest
        if(pdest != MPI_PROC_NULL && iE >= NEloc-iEhbaromegaLO)
        {
          iEplus = iE - (NEloc - iEhbaromegaLO);
          MPI_Isend(&Gless[(iK*NEloc+iEplus)*Np*Np],Np*Np,MPI_DOUBLE_COMPLEX,pdest,42,comm2d,&rqE[2]);
        }
        // recv from source
        if(psource != MPI_PROC_NULL && iE >= NEloc-iEhbaromegaLO)
        {
          MPI_Irecv(rbuff2,Np*Np,MPI_DOUBLE_COMPLEX,psource,42,comm2d,&rqE[3]);
          // pointer
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
      
      
      
      // Update local
      #pragma omp parallel for private(iQ,jz,iz,iQglo,im) collapse(2)
      for( iQ=0; iQ<NKloc; iQ++ )
      {
        for( jz=0; jz<Np; jz++ )
        {
          iQglo=iQ+coords[0]*NKloc;
          for( iz=0; iz<Np; iz++ )
          {
            im = abs(iz-jz);
            sigless[(iQ*NEloc+iE)*Np*Np+jz*Np+iz] += dK * qe * hbaromegaLO/(4.0*pow(M_PI,2.0)) * 
              (1.0/(epsinfr*eps0) - 1.0/(eps0r*eps0)) * K[iKglo] * F[ im+ Np*iKglo + Np*NK*iQglo ] * 
              Mtilde * (fac1 * pbuff1[jz*Np+iz] + fac2 * pbuff2[jz*Np+iz]); 
          }
        }
      }
      
      
      if(dims[0] > 1)
      {
        //printf("Communication for Q integration\n");
        #pragma omp parallel for private(iz,jz)
        for( jz=0; jz<Np; jz++ )
          for( iz=0; iz<Np; iz++ )
            sbuffH[jz*Np+iz] = fac1 * pbuff1[jz*Np+iz] + fac2 * pbuff2[jz*Np+iz];
      
      
        for( in=1; in<dims[0]; in++ )
        {
          MPI_Cart_shift( comm2d_per, 0, in, &hsource, &hdest );
          
          MPI_Cart_coords(comm2d,hsource,ndims,coordsH);
          iKglo2=iK+coordsH[0]*NKloc;
          
          // send to in+1%dims[0] = hdest
          MPI_Isend(sbuffH,Np*Np,MPI_DOUBLE_COMPLEX,hdest,43,comm2d_per,&rqH[0]);
          // recv from in+dims[0]%dims[0] = hsource
          MPI_Irecv(rbuffH,Np*Np,MPI_DOUBLE_COMPLEX,hsource,43,comm2d_per,&rqH[1]);
          
          MPI_Wait(&rqH[1],&statusH[1]);
          
          #pragma omp parallel for private(iQ,jz,iz,iQglo,im) collapse(2)
          for( iQ=0; iQ<NKloc; iQ++ )
          {
            for( jz=0; jz<Np; jz++ )
            {
              iQglo=iQ+coords[0]*NKloc;
              for( iz=0; iz<Np; iz++ )
              {
                im = abs(iz-jz);
                sigless[(iQ*NEloc+iE)*Np*Np+jz*Np+iz] += dK * qe * hbaromegaLO/(4.0*pow(M_PI,2.0)) * 
                  (1.0/(epsinfr*eps0) - 1.0/(eps0r*eps0)) * K[iKglo2] * F[ im+ Np*iKglo2 + Np*NK*iQglo ] * 
                  Mtilde * rbuffH[jz*Np+iz];
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
