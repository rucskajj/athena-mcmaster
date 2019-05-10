#include "../copyright.h"
/*============================================================================*/
/*! \file selfg_fft_par.c
 *  \brief Contains functions to solve Poisson's equation for self-gravity of 
 *   the gas as well as particles in 1D, 2D and 3D using FFTs (actually, the 1D  *   algorithm uses Forward Elimination followed by Back Substitution: FEBS).
 *
 *   These functions require PERIODIC BCs and use the Jeans swindle.
 *
 *   The 2D and 3D f'ns use FFTW3.x, and for MPI parallel use Steve Plimpton's
 *   block decomposition routines added by N. Lemaster to /athena/fftsrc.
 *   This means to use these fns the code must be
 *   - (1) configured with --with-gravity=fft_par --enable-fft
 *   - (2) compiled with links to FFTW libraries (may need to edit Makeoptions)
 *
 *   For NON-PERIODIC BCs, use selfg_multig() functions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - selfg_fft_1d() - actually uses FEBS
 * - selfg_fft_2d() - 2D Poisson solver using FFTs
 * - selfg_fft_3d() - 3D Poisson solver using FFTs
 * - selfg_fft_2d_init() - initializes FFT plans for 2D
 * - selfg_fft_3d_init() - initializes FFT plans for 3D */
/*============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"
#include "../particles/prototypes.h"
#include "../particles/particle.h"

#ifdef SELF_GRAVITY_USING_FFT_PAR

#ifndef FFT_ENABLED
#error self gravity with FFT requires configure --enable-fft
#endif /* FFT_ENABLED */

#ifndef PARTICLES
#error : self gravity with particles requires particles to be enabled.
#endif /* PARTICLES */

/* plans for forward and backward FFTs; work space for FFTW */
static struct ath_2d_fft_plan *fplan2d, *bplan2d;
static struct ath_3d_fft_plan *fplan3d, *bplan3d;
static ath_fft_data *work=NULL;

#ifdef STATIC_MESH_REFINEMENT
#error self gravity with FFT not yet implemented to work with SMR
#endif


/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_3d(DomainS *pD)
 *  \brief Only works for uniform grid, periodic boundary conditions
 */

void selfg_fft_par_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2),dx3sq=(pG->dx3*pG->dx3);
  Real dkx,dky,dkz,pcoeff;

  Real3Vect cell1;              /* one over dx1, dx2, dx3 */
  /* cell1 is a shortcut expressions as well as dimension indicator */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;

  particle_to_grid(pD, property_all);

#ifdef SHEARING_BOX
  Real qomt,Lx,Ly,dt;
  Real kxtdx;
  Real xmin,xmax;
  int ip,jp;
  int nx3=pG->Nx[2]+2*nghost;
  int 
nx2=pG->Nx[1]+2*nghost;
  int nx1=pG->Nx[0]+2*nghost;
  Real ***RollDen, ***UnRollPhi;

  if((RollDen=(Real***)calloc_3d_array(nx3,nx1,nx2,sizeof(Real)))==NULL)
    ath_error("[selfg_fft_3d]: malloc returned a NULL pointer\n");
  if((UnRollPhi=(Real***)calloc_3d_array(nx3,nx1,nx2,sizeof(Real)))==NULL)
    ath_error("[selfg_fft_3d]: malloc returned a NULL pointer\n");

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  dt = pG->time-((int)(qshear*Omega_0*pG->time*Lx/Ly))*Ly/(qshear*Omega_0*Lx);
  qomt = qshear*Omega_0*dt;
#endif

  int il, jl, kl;
  il = pG->is;
  jl = pG->js;
  kl = pG->ks;

/* Copy current potential into old */

  for (k=ks-nghost; k<=ke+nghost; k++){
  for (j=js-nghost; j<=je+nghost; j++){
    for (i=is-nghost; i<=ie+nghost; i++){
      pG->Phi_old[k][j][i] = pG->Phi[k][j][i];
#ifdef SHEARING_BOX
      RollDen[k][i][j] = pG->U[k][j][i].d + pG->Coup[k][j][i].grid_d;
#endif
    }
  }}

/* Forward FFT of 4\piG*(d-d0) */

/* For shearing-box, need to roll density to the nearest periodic point */
#ifdef SHEARING_BOX
  RemapVar(pD,RollDen,-dt);
#endif

  for (k=ks; k<=ke; k++){
  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
#ifdef SHEARING_BOX
        RollDen[k][i][j] - grav_mean_rho;
#else
        pG->U[k][j][i].d - grav_mean_rho + pG->Coup[k][j][i].grid_d; 
#endif
      /*ath_pout(0,"[fft_par] %d %d %d: densities: %g %g %g\n", k,j,i, \
                                    RollDen[k][j][i], \
                                    grav_mean_rho, \
                                    RollDen[k][j][i]-grav_mean_rho);
      */

      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;
    }
  }}

#ifndef SHEARING_BOX 
#ifdef STAR_PARTICLE
   assign_starparticles_3d(pD,work); 
#endif /* STAR_PARTICLE */
#endif /* SHEARING_BOX  */

    
   ath_3d_fft(fplan3d, work);
     
/* Compute potential in Fourier space.  Multiple loops are used to avoid divide
 * by zero at i=is,j=js,k=ks, and to avoid if statement in loop   */
/* To compute kx,ky,kz, note that indices relative to whole Domain are needed */

  dkx = 2.0*PI/(double)(pD->Nx[0]);
  dky = 2.0*PI/(double)(pD->Nx[1]);
  dkz = 2.0*PI/(double)(pD->Nx[2]);

#ifdef SHEARING_BOX
  ip=KCOMP(0,pG->Disp[0],pD->Nx[0]);
  jp=KCOMP(0,pG->Disp[1],pD->Nx[1]);
  kxtdx  = (ip+qomt*Lx/Ly*jp)*dkx;
#endif

  if ((pG->Disp[2])==0 && (pG->Disp[1])==0 && (pG->Disp[0])==0) {
    work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 0.0;
    work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;
  } else {
#ifdef SHEARING_BOX
    pcoeff = 1.0/(((2.0*cos( kxtdx           )-2.0)/dx1sq) +
                  ((2.0*cos((pG->Disp[1])*dky)-2.0)/dx2sq) +
                  ((2.0*cos((pG->Disp[2])*dkz)-2.0)/dx3sq));
#else
    pcoeff = 1.0/(((2.0*cos((pG->Disp[0])*dkx)-2.0)/dx1sq) +
                  ((2.0*cos((pG->Disp[1])*dky)-2.0)/dx2sq) +
                  ((2.0*cos((pG->Disp[2])*dkz)-2.0)/dx3sq));
#endif
    work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
    work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
  }


  for (k=ks+1; k<=ke; k++){
#ifdef SHEARING_BOX
    pcoeff = 1.0/(((2.0*cos( kxtdx                    )-2.0)/dx1sq) +
                  ((2.0*cos((        pG->Disp[1] )*dky)-2.0)/dx2sq) +
                  ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#else
    pcoeff = 1.0/(((2.0*cos((        pG->Disp[0] )*dkx)-2.0)/dx1sq) +
                  ((2.0*cos((        pG->Disp[1] )*dky)-2.0)/dx2sq) +
                  ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#endif
    work[F3DI(0,0,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
    work[F3DI(0,0,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
  }

  for (j=js+1; j<=je; j++){
    for (k=ks; k<=ke; k++){
#ifdef SHEARING_BOX
      jp=KCOMP(j-js ,pG->Disp[1],pD->Nx[1]);
      kxtdx  = (ip+qomt*Lx/Ly*jp)*dkx;
      pcoeff = 1.0/(((2.0*cos( kxtdx                    )-2.0)/dx1sq) +
                    ((2.0*cos(( (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
                    ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#else
      pcoeff = 1.0/(((2.0*cos((        pG->Disp[0] )*dkx)-2.0)/dx1sq) +
                    ((2.0*cos(( (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
                    ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#endif
      work[F3DI(0,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
      work[F3DI(0,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
    }
  }

  for (i=is+1; i<=ie; i++){
  for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
#ifdef SHEARING_BOX
      ip=KCOMP(i-is ,pG->Disp[0],pD->Nx[0]);
      jp=KCOMP(j-js ,pG->Disp[1],pD->Nx[1]);
      kxtdx  = (ip+qomt*Lx/Ly*jp)*dkx;
      pcoeff = 1.0/(((2.0*cos( kxtdx                    )-2.0)/dx1sq) +
                    ((2.0*cos(( (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
                    ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#else
      pcoeff = 1.0/(((2.0*cos(( (i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
                    ((2.0*cos(( (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
                    ((2.0*cos(( (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
#endif
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
    }
  }}

/* Backward FFT and set potential in real space.  Normalization of Phi is over
 * total number of cells in Domain */

  ath_3d_fft(bplan3d, work);

  for (k=ks; k<=ke; k++){
  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
#ifdef SHEARING_BOX
      UnRollPhi[k][i][j] = 
       four_pi_G*work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
        / bplan3d->gcnt;
#else
      pG->Phi[k][j][i] =
       four_pi_G*work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
        / bplan3d->gcnt;
#endif
    }
  }}

#ifdef SHEARING_BOX
  RemapVar(pD,UnRollPhi,dt);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
         pG->Phi[k][j][i] = UnRollPhi[k][i][j];
         pG->GradPhiX1[k][j][i] = -0.5*cell1.x1*(pG->Phi[k][j][i+1]-pG->Phi[k][j][i-1]);
         pG->GradPhiX2[k][j][i] = -0.5*cell1.x2*(pG->Phi[k][j+1][i]-pG->Phi[k][j-1][i]);
         pG->GradPhiX3[k][j][i] = -0.5*cell1.x3*(pG->Phi[k+1][j][i]-pG->Phi[k-1][j][i]);
      }
    }
  }

  free_3d_array(RollDen);
  free_3d_array(UnRollPhi);
#endif

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_3d_init(MeshS *pM)
 *  \brief Initializes plans for forward/backward FFTs, and allocates memory 
 *  needed by FFTW.
 */

void selfg_fft_par_3d_init(MeshS *pM)
{
  DomainS *pD;
  int nl,nd;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        fplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
        bplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
        work = ath_3d_fft_malloc(fplan3d);
      }
    }
  }
}

#endif /* SELF_GRAVITY_USING_FFT_PAR */
