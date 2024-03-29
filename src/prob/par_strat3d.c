#include "copyright.h"
/*============================================================================*/
/*! \file par_strat3d.c
 *  \brief Problem generator for non-linear streaming instability in stratified
 *   disks.
 *
 * PURPOSE: Problem generator for non-linear streaming instability in stratified
 *   disks. This code works in 3D ONLY. Isothermal eos is assumed, and the value
 *   etavk/iso_sound is fixed. MPI domain decomposition in x is allowed, but 
 *   not in z.
 *
 * Perturbation modes:
 * -  ipert = 0: multi-nsh equilibtium
 * -  ipert = 1: white noise within the entire grid
 * -  ipert = 2: non-nsh velocity
 *
 *  Should be configured using --enable-shearing-box and --with-eos=isothermal.
 *  FARGO is recommended.
 *
 * Reference:
 * - Johansen & Youdin, 2007, ApJ, 662, 627
 * - Bai & Stone, 2009, in preparation */
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

#ifndef SHEARING_BOX
#error : The streaming3d problem requires shearing-box to be enabled.
#endif /* SHEARING_BOX */

#ifndef PARTICLES
#error : The streaming3d problem requires particles to be enabled.
#endif /* PARTICLES */

#ifndef ISOTHERMAL
#error : The streaming3d problem requires isothermal equation of state.
#endif /* ISOTHERMAL */

/*------------------------ filewide global variables -------------------------*/
/* flow properties */
Real vsc1,vsc2;
int ipert;
/* domain size variables */
Real x1min,x1max,x2min,x2max,x3min,x3max,Lx,Ly,Lz,Lg;
long Npar;
/* output variables */
long ntrack;   /* number of tracer particles */
long nlis;     /* number of output particles for list output */
int mytype;    /* specific particle type for output */
Real dpar_thresh; /* threshold particle density */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()            - random number generator
 * Normal()          - normal distribution generator
 * Erf()             - error function
 * read_custom_gsd   - read in custom grain size distribution data
 * MultiNSH()        - multiple component NSH equilibrium solver
 * ShearingBoxPot()  - shearing box tidal gravitational potential
 * hst_rho_Vx_dVy()  - total Reynolds stress for history dump
 * property_???()    - particle property selection function
 *============================================================================*/

double ran2(long int *idum);
double Normal(long int *idum);
Real Erf(Real z);

void read_custom_gsd(int n, Real *tstop, Real *weights);
void MultiNSH(int n, Real *tstop, Real *mratio, Real etavk,
                     Real *uxNSH, Real *uyNSH, Real *wxNSH, Real *wyNSH);
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real StratifiedDisk(const Real x1, const Real x2, const Real x3);

static int property_limit(const GrainS *gr, const GrainAux *grsub);
static int property_trace(const GrainS *gr, const GrainAux *grsub);
static int property_type(const GrainS *gr, const GrainAux *grsub);
static int property_dense(const GrainS *gr, const GrainAux *grsub);

extern Real expr_dpar(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V1par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V2par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V3par(const GridS *pG, const int i, const int j, const int k);
extern Real expr_V2(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k,ks,pt,tsmode,gsdmode,weightmode;
  long p,q;
  Real ScaleHg,tsmin,tsmax,tscrit,amin,amax,Hparmin,Hparmax;
  Real *ep,*gsd_weights,*ScaleHpar,epsum,mratio,pwind,rhoaconv,etavk;
  Real *epsilon,*uxNSH,*uyNSH,**wxNSH,**wyNSH;
  Real rhog,h,x1,x2,x3,t,x1p,x2p,x3p,zmin,zmax,dx3_1,b;
  long int iseed, iseedpert; /* Initialize on the first call to ran2 */

  iseedpert = par_geti_def("particle","iseedpert",0);
  iseed = myID_Comm_world+iseedpert;
  /* Initialize on the first call to ran2 */


  if (pDomain->Nx[2] == 1) {
    ath_error("[par_strat3d]: par_strat3d only works for 3D problem.\n");
  }
  
#ifdef MPI_PARALLEL
  if (pDomain->NGrid[2] > 2) {
    ath_error(   
  "[par_strat3d]: The z-domain can not be decomposed into more than 2 grids\n");
  }
#endif

/* Initialize boxsize */
  x1min = pGrid->MinX[0];
  x1max = pGrid->MaxX[0];
  Lx = x1max - x1min;

  x2min = pGrid->MinX[1];
  x2max = pGrid->MaxX[1];
  Ly = x2max - x2min;

  x3min = par_getd("domain1","x3min");
  x3max = par_getd("domain1","x3max");
  Lz = x3max - x3min;

  Lg = nghost*pGrid->dx3; /* size of the ghost zone */

  ks = pGrid->ks;

/* Read initial conditions */
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  ipert = par_geti_def("problem","ipert",1);
  vsc1 = par_getd_def("problem","vsc1",0.05); /* in unit of iso_sound (N.B.!) */
  vsc2 = par_getd_def("problem","vsc2",0.0);

#ifdef SELF_GRAVITY
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = 0.0; //rho0 * (1+mratio);
#endif /* SELF_GRAVITY */

  vsc1 = vsc1 * Iso_csound;
  vsc2 = vsc2 * Iso_csound;

  ScaleHg = Iso_csound/Omega_0;

  /* particle number */
  Npar  = (long)(par_geti("particle","parnumgrid"));

  pGrid->nparticle = Npar*npartypes;
  for (i=0; i<npartypes; i++)
    grproperty[i].num = Npar;

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  ep = (Real*)calloc_1d_array(npartypes, sizeof(Real));
  gsd_weights = (Real*)calloc_1d_array(npartypes, sizeof(Real));
  ScaleHpar = (Real*)calloc_1d_array(npartypes, sizeof(Real));

  epsilon = (Real*)calloc_1d_array(npartypes, sizeof(Real));
  wxNSH   = (Real**)calloc_2d_array(pGrid->Nx[2]+1, npartypes,sizeof(Real));
  wyNSH   = (Real**)calloc_2d_array(pGrid->Nx[2]+1, npartypes,sizeof(Real));
  uxNSH   = (Real*)calloc_1d_array(pGrid->Nx[2]+1, sizeof(Real));
  uyNSH   = (Real*)calloc_1d_array(pGrid->Nx[2]+1, sizeof(Real));

  /* particle stopping time */
  tsmode = par_geti("particle","tsmode");
  gsdmode = par_geti_def("particle","gsdmode", 1);
  weightmode = par_geti_def("particle","weightmode", 1);

  if (tsmode == 3) {/* fixed stopping time */
    tscrit= par_getd("problem","tscrit");
    if (gsdmode == 1) { /* uniform grain size dist */
      tsmin = par_getd("problem","tsmin"); /* in code unit */
      tsmax = par_getd("problem","tsmax");

      for (i=0; i<npartypes; i++) {
        tstop0[i] = tsmin*exp(i*log(tsmax/tsmin)/MAX(npartypes-1,1.0));
        grproperty[i].rad = tstop0[i];
        /* use fully implicit integrator for well coupled particles */
        if (tstop0[i] < tscrit) grproperty[i].integrator = 3;
      }
    }
    if (gsdmode == 2) { /* custom grain size dist */
      read_custom_gsd(npartypes, tstop0, gsd_weights);
  
      for(int i=0; i<npartypes; i++) {
        grproperty[i].rad = tstop0[i];
        /* use fully implicit integrator for well coupled particles */
        if (tstop0[i] < tscrit) grproperty[i].integrator = 3;
      }
    }
  }
  else { 
    if (gsdmode != 1) {
      ath_error("[par_strat3d]: Custom sizes only implemented for tsmode=3.\n");
    }

    amin = par_getd("problem","amin");
    amax = par_getd("problem","amax");

    for (i=0; i<npartypes; i++)
      grproperty[i].rad = amin*exp(i*log(amax/amin)/MAX(npartypes-1,1.0));

    if (tsmode <= 2) {/* Epstein/General regime */
      /* conversion factor for rhoa */
      rhoaconv = par_getd_def("problem","rhoaconv",1.0);

      for (i=0; i<npartypes; i++)
        grrhoa[i]=grproperty[i].rad*rhoaconv;
    }

    if (tsmode == 1)  /* General drag formula */
      alamcoeff = par_getd("problem","alamcoeff");
  }

  if (gsdmode == 2) { /* custom grain size dist */
    /* particle scale height */
    Hparmin = par_getd("problem","hparmin");
    Hparmax = par_getd("problem","hparmax");
    if (Hparmax < Hparmin) {
      ath_error("[par_strat3d]: must have hparmax > hparmin!\n");
    }

    float Hparfac = (Hparmax/Hparmin) - 1.0;
    /* This produces a negative linear relationship between Hpar and tstop */
    for (i=0; i<npartypes; i++)
      ScaleHpar[i] = Hparmin * ((Hparfac+1.0) - 
            Hparfac*(tstop0[i]/tstop0[npartypes-1]));
      //ath_pout(0,"tstop[%d]=%e, ScaleHpar[%d]=%e\n",
      //i, tstop0[i], i, ScaleHpar[i]);
  }
  else {
    /* particle scale height */
    Hparmax = par_getd("problem","hparmax"); /* in unit of gas scale height */
    Hparmin = par_getd("problem","hparmin");
    for (i=0; i<npartypes; i++) 
      ScaleHpar[i] = Hparmax*
                   exp(-i*log(Hparmax/Hparmin)/MAX(npartypes-1,1.0));
  }


#ifdef FEEDBACK
  mratio = par_getd_def("problem","mratio",0.0); /* total mass fraction */
  pwind = par_getd_def("problem","pwind",0.0);   /* power law index */
  if (mratio < 0.0)
    ath_error("[par_strat2d]: mratio must be positive!\n");

  if (gsdmode == 1){ /* uniform grain size dist */
    epsum = 0.0;
    for (i=0; i<npartypes; i++)
    {
      ep[i] = pow(grproperty[i].rad,pwind);	epsum += ep[i];
    }

    for (i=0; i<npartypes; i++)
    {
      ep[i] = mratio*ep[i]/epsum;
      grproperty[i].m = sqrt(2.0*PI)*ScaleHg/Lz*ep[i]*
                                   pGrid->Nx[0]*pGrid->Nx[1]*pGrid->Nx[2]/Npar;
    }
  }
  if (gsdmode == 2 && weightmode == 1){ /* custom grain size dist */
                                        /* using the .m mass parameter */
    epsum = 0.0;
    for (i=0; i<npartypes; i++)
    {
      ep[i] = gsd_weights[i]; epsum += ep[i];
    }

    for (i=0; i<npartypes; i++)
    {
      ep[i] = mratio*ep[i]/epsum;
      grproperty[i].m = sqrt(2.0*PI)*ScaleHg/Lz*ep[i]*
                                   pGrid->Nx[0]*pGrid->Nx[1]*pGrid->Nx[2]/Npar;
    }
  }


#else
 if (gsdmode != 1) {
      ath_error("[par_strat3d]: Custom sizes only implemented for feedback on.\n");
 }
 mratio = 0.0;
  for (i=0; i<npartypes; i++)
    ep[i] = 0.0;
#endif

  /* NSH equilibrium */
  for (k=pGrid->ks; k<=pGrid->ke+1; k++) {

    h = pGrid->MinX[2] + (k-pGrid->ks)*pGrid->dx3;
    q = k - ks;
    etavk = fabs(vsc1+vsc2*SQR(h));

    for (i=0; i<npartypes; i++) {
      epsilon[i] = ep[i]/ScaleHpar[i]*exp(-0.5*SQR(h/ScaleHg)
         *(SQR(1.0/ScaleHpar[i])-1.0))/erf(Lz/(sqrt(8.0)*ScaleHpar[i]*ScaleHg));

      if (tsmode != 3)
        tstop0[i] = get_ts(pGrid,i,exp(-0.5*SQR(h/ScaleHg)),Iso_csound,etavk);
    }

    MultiNSH(npartypes, tstop0, epsilon, etavk,
                              &uxNSH[q], &uyNSH[q], wxNSH[q], wyNSH[q]);
  }

/* Now set initial conditions for the gas */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

    rhog = exp(-0.5*SQR(x3/ScaleHg));
    pGrid->U[k][j][i].d = rhog;

    if (ipert != 1) {/* NSH velocity */
      pGrid->U[k][j][i].M1 = 0.5*rhog*(uxNSH[k-ks]+uxNSH[k-ks+1]);
      pGrid->U[k][j][i].M2 = 0.5*rhog*(uyNSH[k-ks]+uyNSH[k-ks+1]);
    } else {
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
    }

    pGrid->U[k][j][i].M3 = 0.0;
#ifndef FARGO
    pGrid->U[k][j][i].M2 -= qshear*rhog*Omega_0*x1;
#endif

  }}}

/* Now set initial conditions for the particles */
  p = 0;
  dx3_1 = 1.0/pGrid->dx3;
  zmin = pGrid->MinX[2];
  zmax = pGrid->MaxX[2];

  for (pt=0; pt<npartypes; pt++) {

    if ( gsdmode == 2 && weightmode ==2 ){ /* do custom dist via Nparticles */
      Npar *= gsd_weights[pt];
    }

    for (q=0; q<Npar; q++) {

      x1p = x1min + Lx*ran2(&iseed);
      x2p = x2min + Ly*ran2(&iseed);
      x3p = ScaleHpar[pt]*ScaleHg*Normal(&iseed);
      while ((x3p >= zmax) || (x3p < zmin))
        x3p = ScaleHpar[pt]*ScaleHg*Normal(&iseed);

      pGrid->particle[p].property = pt;
      pGrid->particle[p].x1 = x1p;
      pGrid->particle[p].x2 = x2p;
      pGrid->particle[p].x3 = x3p;

      if (ipert != 1) {/* NSH velocity */

        cellk(pGrid, x3p, dx3_1, &k, &b);
        k = k-pGrid->ks;  b = b - pGrid->ks;

        pGrid->particle[p].v1 = (k+1-b)*wxNSH[k][pt]+(b-k)*wxNSH[k+1][pt];
        pGrid->particle[p].v2 = (k+1-b)*wyNSH[k][pt]+(b-k)*wyNSH[k+1][pt];

      } else {

        pGrid->particle[p].v1 = 0.0;
        pGrid->particle[p].v2 = vsc1+vsc2*SQR(x2p);

      }

      pGrid->particle[p].v3 = 0.0;
#ifndef FARGO
      pGrid->particle[p].v2 -= qshear*Omega_0*x1p;
#endif

      pGrid->particle[p].pos = 1; /* grid particle */
      pGrid->particle[p].my_id = p;
#ifdef MPI_PARALLEL
      pGrid->particle[p].init_id = myID_Comm_world;
#endif
      p++;
  }}

/* enroll gravitational potential function, shearing sheet BC functions */
  ShearingBoxPot = StratifiedDisk;

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");

  /* set the # of the particles in list output
   * (by default, output 1 particle per cell)
   */
  nlis = par_geti_def("problem","nlis",pGrid->Nx[0]*pGrid->Nx[1]*pGrid->Nx[2]);

  /* set the number of particles to keep track of */
  ntrack = par_geti_def("problem","ntrack",2000);

  /* set the threshold particle density */
  dpar_thresh = par_geti_def("problem","dpar_thresh",10.0);

  /* finalize */
  free(ep);  free(ScaleHpar);
  free(epsilon);
  free_2d_array(wxNSH);  free_2d_array(wyNSH);
  free(uxNSH);           free(uyNSH);

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * read_custom_gsd         - read in data for custom grain size distribution
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void read_custom_gsd(int n, Real *tstop, Real *weights)
{
  FILE *fp;
  char *gsd_filename = par_gets("particle","cust_gsd_file"); 
  ath_pout(0,"custom_gsd file opening: %s\n", gsd_filename);
  fp = fopen(gsd_filename,"r");
  if (fp == NULL)  
    ath_error("Custom grain size distrubtion file %s could not be opened.\n"
      ,gsd_filename);

  for(int i=0; i<n; i++) {
    fscanf(fp, "%lf", &tstop[i]);
    fscanf(fp, "%lf", &weights[i]);
    //ath_pout(0,"[custom gsd]: %d, %e, %e \n", i, tstop[i], weights[i]);
  }
  return;
}

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  DomainS *pD = (DomainS*)&(pM->Domain[0][0]);
  GridS *pG = pD->Grid;
  ShearingBoxPot = StratifiedDisk;

  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  ipert = par_geti_def("problem","ipert",1);

#ifdef SELF_GRAVITY
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = 0.0; //rho0 * (1+mratio);
#endif /* SELF_GRAVITY */

  x1min = pG->MinX[0];
  x1max = pG->MaxX[0];
  Lx = x1max - x1min;

  x2min = pG->MinX[1];
  x2max = pG->MaxX[1];
  Ly = x2max - x2min;

  x3min = pM->RootMinX[2];
  x3max = pM->RootMaxX[2];
  Lz = x3max - x3min;

  Lg = nghost*pG->dx3; /* size of the ghost zone */

  vsc1 = par_getd_def("problem","vsc1",0.05); /* in unit of iso_sound (N.B.!) */
  vsc2 = par_getd_def("problem","vsc2",0.0);

  vsc1 = vsc1 * Iso_csound;
  vsc2 = vsc2 * Iso_csound;

  Npar  = (int)(sqrt(par_geti("particle","parnumgrid")));
  nlis = par_geti_def("problem","nlis",pG->Nx[0]*pG->Nx[1]*pG->Nx[2]);
  ntrack = par_geti_def("problem","ntrack",2000);

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");

  return;
}

static Real hst_rho_Vx_dVy(const GridS *pG,
                           const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M1*pG->U[k][j][i].M2/pG->U[k][j][i].d;
#else
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d
                            + qshear*Omega_0*x1);
#endif
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  if (strcmp(name,"limit")==0)    return property_limit;
  if (strcmp(name,"trace")==0)    return property_trace;
  if (strcmp(name,"type")==0)     return property_type;
  if (strcmp(name,"dense")==0)    return property_dense;

  return NULL;
}

void gasvshift(const Real x1, const Real x2, const Real x3,
                                    Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Real3Vect *ft, const Real x1, const Real x2, const Real x3,
                                    const Real v1, const Real v2, const Real v3)
{
  Real z,fac;

  ft->x1 -= 2.0*(vsc1 + vsc2*SQR(x3))*Omega_0;

  if(x3 > x3max)
    z = x3-Lz;
  else if (x3 < x3min)
    z = x3+Lz;
  else
    z = x3;

  fac = Lg/(0.5*Lz+Lg-fabs(z));
  ft->x3 -= SQR(Omega_0)*z*(1.0-SQR(fac)*fac); /* 3rd order sharp */
//  ft->x3 -= SQR(Omega_0)*z*(1.0-SQR(fac));  /* 2nd order sharp */

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
/*! \fn static Real StratifiedDisk(const Real x1, const Real x2, const Real x3)
 *  \brief shearing box tidal gravitational potential */
static Real StratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real z,z0,phi=0.0;

#ifndef FARGO
  phi -= qshear*SQR(Omega_0*x1);
#endif
  /* Ensure vertical periodicity in ghost zones */
  if(x3 > x3max)
    z=fabs(x3-Lz);
  else if (x3 < x3min)
    z=fabs(x3+Lz);
  else
    z=fabs(x3);

  phi += 0.5*SQR(Omega_0*z);

  /* smooth the potential at domain edges */
  z0 = 0.5*Lz+Lg;
  phi -= SQR(Omega_0*Lg)*Lg*(2.0*z-z0)/(2.0*SQR(z0-z)); /* 3rd order sharp */

//  phi -= SQR(Omega_0*Lg)*(z/(z0-z)+log(z0/(z0-z))); /* 2nd order sharp */

  return phi;
}

/*! \fn static int property_limit(const GrainS *gr, const GrainAux *grsub)
 *  \brief user defined particle selection function (1: true; 0: false) */
static int property_limit(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->my_id<nlis))
    return 1;
  else
    return 0;
}

/*! \fn static int property_trace(const GrainS *gr, const GrainAux *grsub)
 *  \brief user defined particle selection function (1: true; 0: false) */
static int property_trace(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->my_id<ntrack))
    return 1;
  else
    return 0;
}

/*! \fn static int property_type(const GrainS *gr, const GrainAux *grsub)
 *  \brief user defined particle selection function (1: true; 0: false) */
static int property_type(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (gr->property == mytype))
    return 1;
  else
    return 0;
}

/*! \fn static int property_dense(const GrainS *gr, const GrainAux *grsub)
 *  \brief user defined particle selection function (1: true; 0: false) */
static int property_dense(const GrainS *gr, const GrainAux *grsub)
{
  if ((gr->pos == 1) && (grsub->dpar > dpar_thresh))
    return 1;
  else
    return 0;
}

/*----------------------------------------------------------------------------*/
/*! \fn void MultiNSH(int n, Real *tstop, Real *mratio, Real etavk,
 *                   Real *uxNSH, Real *uyNSH, Real *wxNSH, Real *wyNSH)
 *  \brief Multi-species NSH equilibrium 
 *
 * Input: # of particle types (n), dust stopping time and mass ratio array, and 
 *        drift speed etavk.
 * Output: gas NSH equlibrium velocity (u), and dust NSH equilibrium velocity
 *         array (w).
 */
void MultiNSH(int n, Real *tstop, Real *mratio, Real etavk,
                     Real *uxNSH, Real *uyNSH, Real *wxNSH, Real *wyNSH)
{
  int i,j;
  Real *Lambda1,**Lam1GamP1, **A, **B, **Tmp;

  Lambda1 = (Real*)calloc_1d_array(n, sizeof(Real));     /* Lambda^{-1} */
  Lam1GamP1=(Real**)calloc_2d_array(n, n, sizeof(Real)); /* Lambda1*(1+Gamma) */
  A       = (Real**)calloc_2d_array(n, n, sizeof(Real));
  B       = (Real**)calloc_2d_array(n, n, sizeof(Real));
  Tmp     = (Real**)calloc_2d_array(n, n, sizeof(Real));

  /* definitions */
  for (i=0; i<n; i++){
    for (j=0; j<n; j++)
      Lam1GamP1[i][j] = mratio[j];
    Lam1GamP1[i][i] += 1.0;
    Lambda1[i] = 1.0/(tstop[i]+1.0e-16);
    for (j=0; j<n; j++)
      Lam1GamP1[i][j] *= Lambda1[i];
  }

  /* Calculate A and B */
  MatrixMult(Lam1GamP1, Lam1GamP1, n,n,n, Tmp);
  for (i=0; i<n; i++) Tmp[i][i] += 1.0;
  InverseMatrix(Tmp, n, B);
  for (i=0; i<n; i++)
  for (j=0; j<n; j++)
    B[i][j] *= Lambda1[j];
  MatrixMult(Lam1GamP1, B, n,n,n, A);

  /* Obtain NSH velocities */
  *uxNSH = 0.0;  *uyNSH = 0.0;
  for (i=0; i<n; i++){
    wxNSH[i] = 0.0;
    wyNSH[i] = 0.0;
    for (j=0; j<n; j++){
      wxNSH[i] -= B[i][j];
      wyNSH[i] -= A[i][j];
    }
    wxNSH[i] *= 2.0*etavk;
    wyNSH[i] *= etavk;
    *uxNSH -= mratio[i]*wxNSH[i];
    *uyNSH -= mratio[i]*wyNSH[i];
    wyNSH[i] += etavk;
  }

  free(Lambda1);
  free_2d_array(A);         free_2d_array(B);
  free_2d_array(Lam1GamP1); free_2d_array(Tmp);

  return;
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

/*--------------------------------------------------------------------------- */
/*! \fn double Normal(long int *idum)
 *  \brief Normal distribution random number generator */
double Normal(long int *idum)
{
  double Y,X1,X2;

  X1 = ran2(idum);
  X2 = ran2(idum);

  Y = sqrt(-2.0*log(X1+TINY_NUMBER))*cos(2*PI*X2);

  return Y;
}

/*--------------------------------------------------------------------------- */
/*! \fn Real Erf(Real z)
 *  \brief Error function  */
Real Erf(Real z)
{
  /* arrays of the error function y=erf(x) */
  static double x[101]={
        0.000000e+000,  3.783387e-003,  7.709914e-003,  1.178500e-002,  1.601425e-002,  2.040352e-002,
        2.495885e-002,  2.968653e-002,  3.459307e-002,  3.968525e-002,  4.497008e-002,  5.045486e-002,
        5.614715e-002,  6.205480e-002,  6.818596e-002,  7.454909e-002,  8.115295e-002,  8.800667e-002,
        9.511969e-002,  1.025018e-001,  1.101632e-001,  1.181145e-001,  1.263667e-001,  1.349310e-001,
        1.438193e-001,  1.530440e-001,  1.626176e-001,  1.725534e-001,  1.828652e-001,  1.935671e-001,
        2.046738e-001,  2.162008e-001,  2.281639e-001,  2.405796e-001,  2.534651e-001,  2.668380e-001,
        2.807169e-001,  2.951209e-001,  3.100699e-001,  3.255844e-001,  3.416859e-001,  3.583966e-001,
        3.757395e-001,  3.937386e-001,  4.124186e-001,  4.318054e-001,  4.519256e-001,  4.728071e-001,
        4.944786e-001,  5.169701e-001,  5.403124e-001,  5.645379e-001,  5.896800e-001,  6.157732e-001,
        6.428537e-001,  6.709587e-001,  7.001271e-001,  7.303990e-001,  7.618162e-001,  7.944220e-001,
        8.282614e-001,  8.633812e-001,  8.998296e-001,  9.376570e-001,  9.769156e-001,  1.017659e+000,
        1.059945e+000,  1.103830e+000,  1.149376e+000,  1.196644e+000,  1.245701e+000,  1.296614e+000,
        1.349454e+000,  1.404292e+000,  1.461205e+000,  1.520272e+000,  1.581573e+000,  1.645193e+000,
        1.711221e+000,  1.779746e+000,  1.850864e+000,  1.924673e+000,  2.001274e+000,  2.080774e+000,
        2.163281e+000,  2.248909e+000,  2.337778e+000,  2.430008e+000,  2.525728e+000,  2.625070e+000,
        2.728170e+000,  2.835170e+000,  2.946219e+000,  3.061469e+000,  3.181080e+000,  3.305216e+000,
        3.434048e+000,  3.567755e+000,  3.706521e+000,  3.850536e+000,  4.000000e+000
  };
  static double y[101]={
        0.00000000e+000,  4.26907434e-003,  8.69953340e-003,  1.32973284e-002,  1.80686067e-002,  2.30197153e-002,
        2.81572033e-002,  3.34878242e-002,  3.90185379e-002,  4.47565113e-002,  5.07091186e-002,  5.68839404e-002,
        6.32887618e-002,  6.99315688e-002,  7.68205444e-002,  8.39640613e-002,  9.13706742e-002,  9.90491090e-002,
        1.07008250e-001,  1.15257124e-001,  1.23804883e-001,  1.32660778e-001,  1.41834139e-001,  1.51334337e-001,
        1.61170754e-001,  1.71352743e-001,  1.81889576e-001,  1.92790394e-001,  2.04064148e-001,  2.15719527e-001,
        2.27764884e-001,  2.40208149e-001,  2.53056730e-001,  2.66317410e-001,  2.79996226e-001,  2.94098338e-001,
        3.08627885e-001,  3.23587825e-001,  3.38979770e-001,  3.54803790e-001,  3.71058224e-001,  3.87739454e-001,
        4.04841688e-001,  4.22356710e-001,  4.40273635e-001,  4.58578645e-001,  4.77254725e-001,  4.96281391e-001,
        5.15634428e-001,  5.35285634e-001,  5.55202571e-001,  5.75348359e-001,  5.95681482e-001,  6.16155658e-001,
        6.36719759e-001,  6.57317799e-001,  6.77889021e-001,  6.98368078e-001,  7.18685336e-001,  7.38767318e-001,
        7.58537287e-001,  7.77916009e-001,  7.96822665e-001,  8.15175962e-001,  8.32895397e-001,  8.49902691e-001,
        8.66123358e-001,  8.81488386e-001,  8.95935967e-001,  9.09413237e-001,  9.21877939e-001,  9.33299942e-001,
        9.43662512e-001,  9.52963249e-001,  9.61214608e-001,  9.68443923e-001,  9.74692870e-001,  9.80016358e-001,
        9.84480847e-001,  9.88162149e-001,  9.91142807e-001,  9.93509200e-001,  9.95348535e-001,  9.96745927e-001,
        9.97781755e-001,  9.98529475e-001,  9.99054014e-001,  9.99410828e-001,  9.99645625e-001,  9.99794704e-001,
        9.99885782e-001,  9.99939162e-001,  9.99969080e-001,  9.99985060e-001,  9.99993164e-001,  9.99997050e-001,
        9.99998805e-001,  9.99999548e-001,  9.99999841e-001,  9.99999948e-001,  9.99999985e-001
  };
  double z1, g;
  int il, iu, im;
  z1 = fabs(z);
  /* look for the location of z1 in the x array */
  il=0; iu=100;
  while(iu-il>1) {
    im = (iu+il)/2;
    if (x[im] > z1) iu = im;
    else il = im;
  }
  /* linear interpolation */
  g = (y[iu]*(z1-x[il])+y[il]*(x[iu]-z1))/(x[iu]-x[il]);

  g = MIN(g,1.0);
  if (z<0.0) g = -g;
  return g;
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
