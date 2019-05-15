#include "copyright.h"
/*============================================================================*/
/*! \file  selfgXY.c
 *  \brief  Problem generator for linear shearing self gravitating wave test.
 *
 *  Must be configured using --enable-shearing-box and --with-eos=isothermal.
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

#ifndef SHEARING_BOX
#error : The selfgXY problem requires shearing-box to be enabled.
#endif /* SHEARING_BOX */

#ifndef ISOTHERMAL
#error : The selfgXY problem requires isothermal equation of state.
#endif /* ISOTHERMAL */

/*------------------------ filewide global variables -------------------------*/
Real rho0;
/* eigen vector for a streaming instability mode */
Real Reux,Imux,Reuy,Imuy,Reuz,Imuz,Rerho,Imrho,omg,s;
/* perturbation variables */
Real amp, kx, ky;
int nwaveX, nwaveY;
/* output filename */
char name[50];

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ShearingBoxPot()  - shearing box tidal gravitational potential
 * pert_XY()        - perturbation wave form for linear growth rate test
 *============================================================================*/
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real pert_XY(Real fR, Real fI, Real x, Real y, Real t);
extern Real expr_V3(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k,ip,jp,kp,interp,integrator;
  long p;
  Real x1,x2,x3,t,x1l,x1u,x2l,x2u,x3l,x3u,x1p,x2p,x3p;
  Real x1min,x1max,x2min,x2max,x3min,x3max,Lx,Ly,Lz;
  Real rhog,cs,u1,u2,u3;
  long int iseed = -1; /* Initialize on the first call to ran2 */

  if (pDomain->Nx[2] == 1) {
    ath_error("[selfgXY]: streaming3D only works for 3D problem.\n");
  }

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

/* Read initial conditions */
  rho0 = par_getd_def("problem","rho0", 1.0);
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  amp = par_getd_def("problem","amp",0.0);

#ifdef SELF_GRAVITY
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = rho0;
#endif /* SELF_GRAVITY */

  amp = amp;

  Reux  =  par_getd_def("problem","Reux",0.0)*amp;
  Imux  =  par_getd_def("problem","Imux",0.0)*amp;
  Reuy  =  par_getd_def("problem","Reuy",0.0)*amp;
  Imuy  =  par_getd_def("problem","Imuy",0.0)*amp;
  Reuz  =  0.0;
  Imuz  =  0.0;

  Rerho =  par_getd_def("problem","Rerho",0.0)*amp;
  Imrho =  par_getd_def("problem","Imrho",0.0)*amp;

  omg  =  par_getd_def("problem","omg",0.0)*Omega_0;
  s    =  par_getd_def("problem","s"  ,0.0)*Omega_0;
 
  nwaveX = (int)par_getd_def("problem","nwaveX",1);
  nwaveY = (int)par_getd_def("problem","nwaveY",1);

  kx    = 2.0*PI/Lx*nwaveX;
  ky    = 2.0*PI/Ly*nwaveY;
  /* Reset isothermal sound speed */
  cs = 1.0; // set cs to 1.0, box size set in python to vary k*eta*r
  Iso_csound = cs;
  Iso_csound2 = SQR(Iso_csound);

/* Now set initial conditions for the gas */
  t = 0.0;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

    rhog = rho0 * (1.0+pert_XY(Rerho,Imrho,x1,x2,t));
    u1 = pert_XY(Reux,Imux,x1,x2,t);
    u2 = pert_XY(Reuy,Imuy,x1,x2,t);
    u3 = pert_XY(Reuz,Imuz,x1,x2,t);
   
    //rhog = rho0;  u1 = u2 = u3 = 0.0;

    pGrid->U[k][j][i].d = rhog;

    pGrid->U[k][j][i].M1 = rhog * u1;
    pGrid->U[k][j][i].M2 = rhog * u2;

    pGrid->U[k][j][i].M3 = rhog * u3;
#ifndef FARGO
    pGrid->U[k][j][i].M2 -= qshear*rhog*Omega_0*x1;
#endif

  }}}

/* enroll gravitational potential function, shearing sheet BC functions */
  ShearingBoxPot = UnstratifiedDisk;

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
  fwrite(&Reux, sizeof(Real),1,fp); fwrite(&Imux, sizeof(Real),1,fp);
  fwrite(&Reuy, sizeof(Real),1,fp); fwrite(&Imuy, sizeof(Real),1,fp);
  fwrite(&Reuz, sizeof(Real),1,fp); fwrite(&Imuz, sizeof(Real),1,fp);
  fwrite(&Rerho, sizeof(Real),1,fp);fwrite(&Imrho, sizeof(Real),1,fp);
  fwrite(&omg, sizeof(Real),1,fp);  fwrite(&s, sizeof(Real),1,fp);
  fwrite(&Iso_csound, sizeof(Real),1,fp);
  fwrite(&kx, sizeof(Real),1,fp);   fwrite(&ky, sizeof(Real),1,fp);
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  ShearingBoxPot = UnstratifiedDisk;

  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  amp = par_getd_def("problem","amp",0.0);

  fread(&rho0, sizeof(Real),1,fp);  
 
  fread(&Reux, sizeof(Real),1,fp); fread(&Imux, sizeof(Real),1,fp);
  fread(&Reuy, sizeof(Real),1,fp); fread(&Imuy, sizeof(Real),1,fp);
  fread(&Reuz, sizeof(Real),1,fp); fread(&Imuz, sizeof(Real),1,fp);
  fread(&Rerho, sizeof(Real),1,fp);fread(&Imrho, sizeof(Real),1,fp);
  fread(&omg, sizeof(Real),1,fp);  fread(&s, sizeof(Real),1,fp);
  fread(&Iso_csound, sizeof(Real),1,fp);  Iso_csound2 = SQR(Iso_csound);
  fread(&kx, sizeof(Real),1,fp);   fread(&ky, sizeof(Real),1,fp);
  return;
}

/*! \fn static Real expr_rhodif(const GridS *pG, const int i, const int j, 
 *			        const int k)
 *  \brief difd */
static Real expr_rhodif(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return (pG->U[k][j][i].d - rho0)/amp;
}

/*! \fn static Real expr_dVx(const GridS *pG, const int i, const int j, 
 *			     const int k)
 *  \brief dVx */
static Real expr_dVx(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->U[k][j][i].M1/pG->U[k][j][i].d;
}

/*! \fn static Real expr_dVy(const GridS *pG, const int i, const int j, 
 *			     const int k)
 *  \brief dVy */
static Real expr_dVy(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M2/pG->U[k][j][i].d;
#else
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
}

ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"difd")==0) return expr_rhodif;
  if(strcmp(expr,"dVx")==0) return expr_dVx;
  if(strcmp(expr,"dVy")==0) return expr_dVy;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  return NULL;
}

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
/*! \fn static Real UnstratifiedDisk(const Real x1, const Real x2,const Real x3)
 *  \brief shearing box tidal gravitational potential */
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*SQR(Omega_0*x1);
#endif
  return phi;
}

/*! \fn static Real pert_XY(Real fR, Real fI, Real x, Real y, Real t)
 *  \brief shearing planar perturbation mode */
static Real pert_XY(Real fR, Real fI, Real x, Real y, Real t)
{
  return (fR*cos(kx*x-ky*y-omg*t)-fI*sin(kx*x-ky*y-omg*t));
}


