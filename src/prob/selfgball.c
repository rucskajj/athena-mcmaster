#include "copyright.h"
/*============================================================================*/
/*! \file  selfgball.c
 *  \brief  Problem generator for self gravitating particle ball test.
 */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

#ifndef PARTICLES
#error : The selfgXY problem requires particles to be enabled.
#endif /* PARTICLES */

/*------------------------ filewide global variables -------------------------*/
/* NSH equilibrium parameters */
Real rho0, mratio, etavk;
/* particle number variables */
int Npar,Npar3,downsamp,Nx;
/* perturbation variables */
Real amp, kx, ky;
int ipert, nwaveX, nwaveY;
/* output filename */
char name[50];

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()            - random number generator
 * OutputModeAmplitude() - output the perturbation amplitude
 * ShearingBoxPot()  - shearing box tidal gravitational potential
 * pert_???()        - perturbation wave form for linear growth rate test
 * property_mybin()  - particle property selection function
 *============================================================================*/
double ran2(long int *idum);
void OutputModeAmplitude(MeshS *pM, OutputS *pOut);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static int property_mybin(const GrainS *gr, const GrainAux *grsub);
extern Real expr_V3(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V1par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V2par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V3par(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k,ip,jp,kp,interp,integrator;
  long p;
  Real x1,x2,x3,t,x1l,x1u,x2l,x2u,x3l,x3u,x1p,x2p,x3p,rp;
  Real x1min,x1max,x2min,x2max,x3min,x3max,Lx,Ly,Lz;
  Real rhog,parrad,cs;
  long int iseed = -1; /* Initialize on the first call to ran2 */

  if (pDomain->Nx[2] == 1) {
    ath_error("[selfgball]: selfgball only works for 3D problem.\n");
  }

  interp = par_geti("particle","interp");


/* Initialize boxsize */
  x1min = pDomain->RootMinX[0];
  x1max = pDomain->RootMaxX[0];
  Lx = x1max - x1min;

  x2min = pDomain->RootMinX[1];
  x2max = pDomain->RootMaxX[1];
  Ly = x2max - x2min;

  x3min = pDomain->RootMinX[2];
  x3max = pDomain->RootMaxX[2];
  Lz = x3max - x3min;

  Nx = pDomain->Nx[0]; /* used for particle selection function */

/* Read initial conditions */
  rho0 = par_getd_def("problem","rho0", 1.0);
  parrad = par_getd_def("problem","parrad", 0.5); /* in units of x1 box size */

  /* particle number */
  if (npartypes != 1)
    ath_error("[selfgball]: This test only allows ONE particle species!\n");

  Npar = 1; Npar3 = 1;

  /* set the down sampling of the particle list output
   * (by default, output 1 particle per cell)
   */
  downsamp = par_geti_def("problem","downsamp",Npar3);

  /* particle stopping time */
  tstop0[0] = par_getd("problem","tstop"); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[selfgball]: This test works only for fixed stopping time!\n");

  /* assign particle effective mass */
  mratio = par_getd_def("problem","mratio",0.0); /* total mass fraction */

#ifdef SELF_GRAVITY
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = 0.0;//rho0 * (1+mratio);
#endif /* SELF_GRAVITY */

#ifdef FEEDBACK
  grproperty[0].m = rho0*mratio/Npar3;
#else
  grav_mean_rho = 0.0;//1.0; // dust particle mass set to 1 in 
                             // particle_to_grid routine with FB off
  mratio = 0.0;
#endif

  integrator = par_geti("particle","integrator");
  if(integrator != 2){
     ath_error("[selfgball]: Currently testing with drag on dust off!\
 integrator must be 2 for the semimp integrator in \
 integrators_particle.c .\n");
  }

  /* Reset isothermal sound speed */
  cs = 1.0; // set cs to 1.0, box size set in python to vary k*eta*r
  Iso_csound = cs;
  Iso_csound2 = SQR(Iso_csound);

/* Now set initial conditions for the gas */
  t = 0.0;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    pGrid->U[k][j][i].d = rho0;

    pGrid->U[k][j][i].M1 = 0.0;
    pGrid->U[k][j][i].M2 = 0.0;
    pGrid->U[k][j][i].M3 = 0.0;

  }}}

  ath_pout(0, "[selfgball] got here 1. \n");


/* Now set initial conditions for the particles */
  p = 0;

  for (k=pGrid->ks; k<=pGrid->ke; k++)
  {
    x3l = pGrid->MinX[2] + (k-pGrid->ks)*pGrid->dx3;
    x3u = pGrid->MinX[2] + (k-pGrid->ks+1.0)*pGrid->dx3;

    for (j=pGrid->js; j<=pGrid->je; j++)
    {
      x2l = pGrid->MinX[1] + (j-pGrid->js)*pGrid->dx2;
      x2u = pGrid->MinX[1] + (j-pGrid->js+1.0)*pGrid->dx2;

      for (i=pGrid->is; i<=pGrid->ie; i++)
      {
        x1l = pGrid->MinX[0] + (i-pGrid->is)*pGrid->dx1;
        x1u = pGrid->MinX[0] + (i-pGrid->is+1.0)*pGrid->dx1;

        for (ip=0;ip<Npar;ip++)
        {
          x1p = x1l+pGrid->dx1*(ip+0.5);

          for (jp=0;jp<Npar;jp++)
          {
            x2p = x2l+pGrid->dx2*(jp+0.5);

            for (kp=0;kp<Npar;kp++)
            {
              x3p = x3l+pGrid->dx3*(kp+0.5);

              rp = sqrt(SQR(x1p) + SQR(x2p) + SQR(x3p));


              if(rp < parrad*(0.5*Lx)){

              //ath_pout(0, "[streaming3d]: rp: %g, %g, %g, %g %g\n",
              //                               x1p, x2p, x3p, rp, parrad*0.5*Lx);

              pGrid->particle[p].property = 0;
              pGrid->particle[p].x1 = x1p;
              pGrid->particle[p].x2 = x2p;
              pGrid->particle[p].x3 = x3p;

              pGrid->particle[p].v1 = 0.0;
              pGrid->particle[p].v2 = 0.0;
              pGrid->particle[p].v3 = 0.0;

              pGrid->particle[p].pos = 1; /* grid particle */
              pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
              pGrid->particle[p].init_id = myID_Comm_world;
#endif
              p += 1;
              }
            }
          }
        }
      }
    }
  }

  pGrid->nparticle = p;
  grproperty[0].num = pGrid->nparticle;

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  if (myID_Comm_world == 0) {
    /* flush output file */
    sprintf(name, "%s_%d_%d.dat", "SelfgXY",Nx,ipert);
    FILE *fid = fopen(name,"w");
    fclose(fid);
#ifdef MPI_PARALLEL
    sprintf(name, "../%s_%d_%d.dat", "SelfgXY",Nx,ipert);
#endif
  }

  ath_pout(0, "[selfgball] got here 2. %d \n", p);

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  fwrite(&rho0, sizeof(Real),1,fp);  fwrite(&mratio, sizeof(Real),1,fp);
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  downsamp = par_geti_def("problem","downsamp",Npar3);

  fread(&rho0, sizeof(Real),1,fp);  fread(&mratio, sizeof(Real),1,fp);

  if (myID_Comm_world == 0)
#ifdef MPI_PARALLEL
    sprintf(name, "../%s_%d_%d.dat", "SelfgXY",Nx,ipert);
#else
    sprintf(name, "%s_%d_%d.dat", "SelfgXY",Nx,ipert);
#endif

  return;
}

/*! \fn static Real expr_GradPhiX1(const GridS *pG, const int i, const int j, 
 *			        const int k)
 *  \brief GradPhiX1 */
static Real expr_GradPhiX1(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->GradPhiX1[k][j][i];
}

/*! \fn static Real expr_GradPhiX2(const GridS *pG, const int i, const int j, 
 *			        const int k)
 *  \brief GradPhiX2 */
static Real expr_GradPhiX2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->GradPhiX2[k][j][i];
}

/*! \fn static Real expr_GradPhiX3(const GridS *pG, const int i, const int j, 
 *			        const int k)
 *  \brief GradPhiX3 */
static Real expr_GradPhiX3(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->GradPhiX3[k][j][i];
}

ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"gphi1")==0) return expr_GradPhiX1;
  if(strcmp(expr,"gphi2")==0) return expr_GradPhiX2;
  if(strcmp(expr,"gphi3")==0) return expr_GradPhiX3;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  return NULL;
}

#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  if (strcmp(name,"downsamp")==0) return property_mybin;
  return NULL;
}

/*! \fn void gasvshift(const Real x1, const Real x2, const Real x3,
 *                                  Real *u1, Real *u2, Real *u3)
 *  \brief Gas velocity shift*/
void gasvshift(const Real x1, const Real x2, const Real x3,
                                    Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Real3Vect *ft, const Real x1, const Real x2, const Real x3,
                                    const Real v1, const Real v2, const Real v3)
{
  //ft->x1 -= 2.0*etavk*Omega_0;
  return;
}
#endif

void Userwork_in_loop(MeshS *pM)
{
  return;
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(MeshS *pM)
{
  return;
}
 

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */

/*! \fn static int property_mybin(const Grain *gr, const GrainAux *grsub)
 *  \brief user defined particle selection function (1: true; 0: false) */
static int property_mybin(const GrainS *gr, const GrainAux *grsub)
{
  long a,b,c,d,e,ds,sp;

  sp = MAX(downsamp/Npar3,1);  /* spacing in cells */
  ds = Npar3*sp;               /* actual dowmsampling */

  a = gr->my_id/ds;
  b = gr->my_id - a*ds;

  c = gr->my_id/(Npar3*Nx);    /* column number */
  d = c/sp;
  e = c-sp*d;

  if ((e == 0) && (b == 0) && (gr->pos == 1))
    return 1;
  else
    return 0;
}


/*------------------------------------------------------------------------------
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  \brief Extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1.
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX
