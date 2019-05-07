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

  //particle_to_grid_fft(pD, property_all);
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

  //ath_pout(0,"[fft_par] indices: [%d][%d][%d] \n", kl,jl,il);
  ath_pout(0,"[fft_par] indices: [%d][%d][%d] \n", ke,je,ie);

  //ath_pout(0,"[partogrid_fft] [%d][%d][%d] 1:  %f \n", kt,jt,it, pG->Coup[kt][jt][it].grid_d);

/* Copy current potential into old */

  for (k=ks-nghost; k<=ke+nghost; k++){
  for (j=js-nghost; j<=je+nghost; j++){
    for (i=is-nghost; i<=ie+nghost; i++){
      pG->Phi_old[k][j][i] = pG->Phi[k][j][i];// pG->Coup[k][j][i].grid_d;
#ifdef SHEARING_BOX
      RollDen[k][i][j] = pG->U[k][j][i].d + pG->Coup[k][j][i].grid_d;

    //  ath_pout(0,"[fft_par] indices: [%d][%d][%d] \n", kl,jl,il);

      //if(abs(pG->Coup[k][j][i].grid_d) < 1e-6){
      //ath_pout(0,"[fft_par] %d %d %d: densities: %g %g\n", k,j,i,pG->U[k][j][i].d, 
      //                              pG->Coup[k][j][i].grid_d);
      //}
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
        RollDen[k][i][j] - grav_mean_rho;// + pG->Coup[k][j][i].grid_d;
#else
        pG->U[k][j][i].d - grav_mean_rho + pG->Coup[k][j][i].grid_d; 
#endif
      /*ath_pout(0,"[fft_par] %d %d %d: densities: %g %g %g\n", k,j,i, \
                                    RollDen[k][j][i], \
                                    grav_mean_rho, \
                                    RollDen[k][j][i]-grav_mean_rho);
      */
7
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

/*----------------------------------------------------------------------------*/
/*! \fn void particle_to_grid(Grid *pG, PropFun_t par_prop)
 *  \brief Bin the particles to grid cells
 */
void particle_to_grid_fft(DomainS *pD, PropFun_t par_prop)
{
  GridS *pG = pD->Grid;
  int i,j,k, is,js,ks, i0,j0,k0, i1,j1,k1, i2,j2,k2;

  /* left and right limit of grid indices */
  //int ilp,iup, jlp,jup, klp,kup;
  /* number of neighbouring cells involved in 1D interpolation */
  int ncell, interp;
  int npartypes;               /*!< number of particle types */
  //Grain_Property *grproperty;  /*!< array of particle properties of all types */

  /* allocate memory for particle properties */
  //grproperty = (Grain_Property*)calloc_1d_array(npartypes,
  //                                            sizeof(Grain_Property));
  //if (grproperty == NULL) ath_error("[init_particle]: Error allocating memory.\n");

  /* set the interpolation function pointer */
  interp = par_geti_def("particle","interp",2);
  if (interp == 1)
  { /* linear interpolation */
    getweight = getwei_linear;
    ncell = 2;
  }
  else if (interp == 2)
  { /* TSC interpolation */
    getweight = getwei_TSC;
    ncell = 3;
  }
  else if (interp == 3)
  { /* Quadratic polynomial interpolation */
    getweight = getwei_QP;
    ncell = 3;
  }
  else
    ath_error("[init_particle]: Value of interp must be 1, 2 or 3!\n");

  int n0 = ncell-1;
  long p;
  Real drho;
  Real weight[3][3][3];
  Real3Vect cell1;
  GrainS *gr;

  /* Get grid limit related quantities */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;
  else                cell1.x1 = 0.0;

  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;
  else                cell1.x2 = 0.0;

  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;
  else                cell1.x3 = 0.0;



  int m1, m2, m3;   /* dimension flags */
  if (pG->Nx[0] > 1) m1 = 1;
  else m1 = 0;

  if (pG->Nx[1] > 1) m2 = 1;
  else m2 = 0;

  if (pG->Nx[2] > 1) m3 = 1;
  else m3 = 0;

/* set left and right grid indices */
/*  ilp = pG->is - m1*nghost;
  iup = pG->ie + m1*nghost;

  jlp = pG->js - m2*nghost;
  jup = pG->je + m2*nghost;

  klp = pG->ks - m3*nghost;
  kup = pG->ke + m3*nghost; */

  //ath_pout(0,"[par2grid_fft] %d %d %d \n", kup, jup, iup);

  /* initialization */
  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++) {

        //ath_pout(0,"[par2grid_fft] Made it in this loop.\n");

        pG->Coup[k][j][i].grid_d = 0.0;
        pG->Coup[k][j][i].grid_v1 = 0.0;
        pG->Coup[k][j][i].grid_v2 = 0.0;
        pG->Coup[k][j][i].grid_v3 = 0.0;

        //ath_pout(0,"[par2grid_fft] %d %f\n", pG->Coup[k][j][i].grid_d);
      }

  //ath_pout(0,"[par2grid_fft] ncell: %d\n", ncell);
  //ath_pout(0,"[par2grid_fft] %d\n", pG->nparticle);

  //ath_pout(2,"[par2grid_fft] %d %g\n", pG->nparticle, pG->Coup[k][j][i].grid_d);
  //ath_pout(3,"[par2grid_fft] %d %g\n", pG->nparticle, pG->Coup[k][j][i].grid_d);


  /* bin the particles */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);

    //ath_pout(0,"[par2grid_fft] In this loop.\n");

    /* judge if the particle should be selected */
    if ((*par_prop)(gr, &(pG->parsub[p]))) {/* 1: true; 0: false */


      //ath_pout(0,"[par2grid_fft] 2: In this loop.\n");

      getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &is, &js, &ks);

      /* distribute particles */
      k1 = MAX(ks, klp);    k2 = MIN(ks+n0, kup);
      j1 = MAX(js, jlp);    j2 = MIN(js+n0, jup);
      i1 = MAX(is, ilp);    i2 = MIN(is+n0, iup);

      ath_pout(0,"[par2grid_fft] nghost: %d \n", nghost);

      ath_pout(0,"[par2grid_fft] up: %d %d %d \n", kup, jup, iup);
      ath_pout(0,"[par2grid_fft] lp: %d %d %d \n", klp, jlp, ilp);

      ath_pout(0,"[par2grid_fft] s:  %d %d %d \n", ks, js, is);
      ath_pout(0,"[par2grid_fft] 1:  %d %d %d \n", k1, j1, i1);
      ath_pout(0,"[par2grid_fft] 2:  %d %d %d \n", k2, j2, i2);  

      for (k=k1; k<=k2; k++) {
        k0 = k-k1;
        for (j=j1; j<=j2; j++) {
          j0 = j-j1;
          for (i=i1; i<=i2; i++) {
            i0 = i-i1;

            //ath_pout(0,"[par2grid_fft] ind:   %d %d %d \n", k, j, i);
            //ath_pout(0,"[par2grid_fft] ind0:  %d %d %d \n", k0, j0, i0);

            //ath_pout(0,"[par2grid_fft] 3: In this loop.\n");

            /* interpolate the particles to the grid */
#ifdef FEEDBACK
            drho = grproperty[gr->property].m;
#else
            drho = 1.0;
#endif
            //ath_pout(0,"[fft_par] %d %d %d: den calc: %g %g\n", k,j,i, weight[k0][j0][i0], drho);
        
            pG->Coup[k][j][i].grid_d  += weight[k0][j0][i0]*drho;
            pG->Coup[k][j][i].grid_v1 += weight[k0][j0][i0]*drho*gr->v1;
            pG->Coup[k][j][i].grid_v2 += weight[k0][j0][i0]*drho*gr->v2;
            pG->Coup[k][j][i].grid_v3 += weight[k0][j0][i0]*drho*gr->v3;
            

            //ath_pout(0,"[par2grid_fft] %g %g\n", drho, pG->Coup[k][j][i].grid_d);

          }
        }
      }
    }
  }

/* deposit ghost zone values into the boundary zones */
  exchange_gpcouple(pD, 0);

  return;
}

#endif /* SELF_GRAVITY_USING_FFT_PAR */
