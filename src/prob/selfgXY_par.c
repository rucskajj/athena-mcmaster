#include "copyright.h"
/*============================================================================*/
/*! \file  selfgXY.c
 *  \brief  Problem generator for linear shearing self gravitating wave test.
 *
 * Perturbation modes:
 *  - ipert = 1: perturbation to gas onky
 *  - ipert = 2: perturbation to dust only
 *  - ipert = 3: pertutbation to gas and dust
 *  - ipert = 4: random perturbation (warm start)
 *
 *  Must be configured using --enable-shearing-box and --with-eos=isothermal.
 *  FARGO is need to establish the NSH equilibrium (ipert=0,1,2).             */
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
#error : The selfgXY problem requires shearing-box to be enabled.
#endif /* SHEARING_BOX */

#ifndef PARTICLES
#error : The selfgXY problem requires particles to be enabled.
#endif /* PARTICLES */

#ifndef ISOTHERMAL
#error : The selfgXY problem requires isothermal equation of state.
#endif /* ISOTHERMAL */

/*------------------------ filewide global variables -------------------------*/
/* NSH equilibrium parameters */
Real rho0, mratio, etavk, uxNSH, uyNSH, wxNSH, wyNSH;
/* particle number variables */
int Npar,Npar3,downsamp,Nx;
/* eigen vector for a streaming instability mode */
Real Reux,Imux,Reuy,Imuy,Reuz,Imuz,Rewx,Imwx,Rewy,Imwy,Rewz,Imwz,Rerho,Imrho,omg,s;
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
static Real pert_even(Real fR, Real fI, Real x, Real z, Real t);
static Real pert_odd(Real fR, Real fI, Real x, Real z, Real t);
static Real pert_XY(Real fR, Real fI, Real x, Real y, Real t);
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
  Real x1,x2,x3,t,x1l,x1u,x2l,x2u,x3l,x3u,x1p,x2p,x3p,paramp,factor2,reduct;
  Real x1min,x1max,x2min,x2max,x3min,x3max,Lx,Ly,Lz;
  Real rhog,cs,u1,u2,u3,w1,w2,w3,denorm1;
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

  Nx = pDomain->Nx[0]; /* used for particle selection function */

/* Read initial conditions */
  rho0 = par_getd_def("problem","rho0", 1.0);
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  amp = par_getd_def("problem","amp",0.0);
  ipert = par_geti_def("problem","ipert", 1);
  etavk = par_getd_def("problem","etavk",0.05);/* in unit of iso_sound (N.B.) */

  /* particle number */
  if (npartypes != 1)
    ath_error("[selfgXY]: This test only allows ONE particle species!\n");

  Npar  = (int)(pow(par_geti("particle","parnumcell"),1.0/3.0));
  Npar3 = Npar*SQR(Npar);

  pGrid->nparticle         = Npar3*pGrid->Nx[0]*pGrid->Nx[1]*pGrid->Nx[2];
  grproperty[0].num = pGrid->nparticle;

  if (pGrid->nparticle+2 > pGrid->arrsize)
    particle_realloc(pGrid, pGrid->nparticle+2);

  /* set the down sampling of the particle list output
   * (by default, output 1 particle per cell)
   */
  downsamp = par_geti_def("problem","downsamp",Npar3);

  /* particle stopping time */
  tstop0[0] = par_getd("problem","tstop"); /* in code unit */
  if (par_geti("particle","tsmode") != 3)
    ath_error("[selfgXY]: This test works only for fixed stopping time!\n");

  /* assign particle effective mass */
  mratio = par_getd_def("problem","mratio",0.0); /* total mass fraction */

#ifdef SELF_GRAVITY
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = rho0 * (1+mratio);
#endif /* SELF_GRAVITY */

#ifdef FEEDBACK
  grproperty[0].m = rho0*mratio/Npar3;
#else
  grav_mean_rho = 1.0; // dust particle mass set to 1 in 
                             // particle_to_grid routine with FB off
  mratio = 0.0;
#endif

  integrator = par_geti("particle","integrator");
  if(integrator != 2){
     ath_error("[selfgXY]: Currently testing with drag on dust off!\
 integrator must be 2 for the semimp integrator in \
 integrators_particle.c .\n");
  }

  Reux  =  par_getd_def("problem","Reux",0.0)*amp;
  Imux  =  par_getd_def("problem","Imux",0.0)*amp;
  Reuy  =  par_getd_def("problem","Reuy",0.0)*amp;
  Imuy  =  par_getd_def("problem","Imuy",0.0)*amp;
  Reuz  =  0.0;
  Imuz  =  0.0;

  Rerho =  par_getd_def("problem","Rerho",0.0)*amp;
  Imrho =  par_getd_def("problem","Imrho",0.0)*amp;

  Rewx  =  par_getd_def("problem","Rewx",0.0)*amp;
  Imwx  =  par_getd_def("problem","Imwx",0.0)*amp;
  Rewy  =  par_getd_def("problem","Rewy",0.0)*amp;
  Imwy  =  par_getd_def("problem","Imwy",0.0)*amp;
  Rewz  =  0.0;
  Imwz  =  0.0;

  omg  =  par_getd_def("problem","omg",0.0)*Omega_0;
  s    =  par_getd_def("problem","s"  ,0.0)*Omega_0;

  /* Adjust code units */ 
  nwaveX = (int)par_getd_def("problem","nwaveX",1);
  nwaveY = (int)par_getd_def("problem","nwaveY",1);

  kx    = 2.0*PI/Lx*nwaveX;
  ky    = 2.0*PI/Ly*nwaveY;

  ath_pout(0, "ky: %f\n", ky);

  /* Reset isothermal sound speed */
  cs = 1.0; // set cs to 1.0, box size set in python to vary k*eta*r
  Iso_csound = cs;
  Iso_csound2 = SQR(Iso_csound);

  interp = par_geti("particle","interp");
  if (interp == 3) {/* QP */
    paramp = amp*kx*pGrid->dx1/sin(kx*pGrid->dx1);
    factor2 = 0.5*tan(kx*pGrid->dx1)/(kx*pGrid->dx1);
    reduct = 1.0;
  }
  else if (interp == 2) {/* TSC */
      paramp = 4.0*amp*kx*pGrid->dx1/sin(kx*pGrid->dx1)/
                                  (3.0+cos(ky*pGrid->dx2));
      factor2 = 0.5*tan(kx*pGrid->dx1)/(kx*pGrid->dx1);
      reduct = 1.0/(1.0-0.25*SQR(kx*pGrid->dx1)); reduct=1.0;
  }
  else {
      paramp = amp;
      factor2 = 0.5;
      reduct = 1.0;
  }
  etavk = etavk * Iso_csound; /* switch to code unit (N.B.!) */

  /* calculate NSH equilibrium velocity */
  denorm1 = 1.0/(SQR(1.0+mratio)+SQR(tstop0[0]*Omega_0));

  wxNSH = -2.0*tstop0[0]*Omega_0*denorm1*etavk;
  wyNSH = -(1.0+mratio)*denorm1*etavk;

  uxNSH = -mratio*wxNSH;
  uyNSH = -mratio*wyNSH;

  wyNSH += etavk;


  /* NOTE. Turning off background velocity stuff. */
  uxNSH = 0.0; uyNSH = 0.0; wxNSH = 0.0, wyNSH = 0.0;
  etavk = 1.0;

  ath_pout(0,"etavk=%f, Iso_csound=%f\n",etavk,Iso_csound);

/* Now set initial conditions for the gas */
  t = 0.0;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

    if ((ipert == 1) || (ipert == 3)) {
      rhog = rho0 * (1.0+reduct*pert_XY(Rerho,Imrho,x1,x2,t));
      u1 = etavk * reduct * pert_XY(Reux,Imux,x1,x2,t);
      u3 = etavk * reduct * pert_XY(Reuz,Imuz,x1,x2,t);
      u2 = etavk * reduct * pert_XY(Reuy,Imuy,x1,x2,t);
    } else {
      rhog = rho0;  u1 = u2 = u3 = 0.0;
    }

    //ath_pout(0,"u1=%g, u2=%g, u3=%g\n",u1,u2,u3);

    pGrid->U[k][j][i].d = rhog;

    pGrid->U[k][j][i].M1 = rhog * (uxNSH+u1);
    pGrid->U[k][j][i].M2 = rhog * (uyNSH+u2);

    pGrid->U[k][j][i].M3 = rhog * u3;
#ifndef FARGO
    pGrid->U[k][j][i].M2 -= qshear*rhog*Omega_0*x1;
#endif

  }}}

/* Now set initial conditions for the particles */
  p = 0;
  Lx = pGrid->Nx[0]*pGrid->dx1;
  x1min = pGrid->MinX[0];

  Ly = pGrid->Nx[1]*pGrid->dx2;
  x2min = pGrid->MinX[1];

  Lz = pGrid->Nx[2]*pGrid->dx3;
  x3min = pGrid->MinX[2];

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

          if (ipert != 4) /* quasi-uniform distribution */
            x1p = x1l+pGrid->dx1/Npar*(ip+0.5);

          for (jp=0;jp<Npar;jp++)
          {
            if (ipert != 4) /* quasi-uniform distribution */
              x2p = x2l+pGrid->dx2/Npar*(jp+0.5);

            for (kp=0;kp<Npar;kp++)
            {
              if (ipert == 4){ /* ramdom particle position in the grid */
                x1p = x1min + Lx*ran2(&iseed);
                x2p = x2min + Ly*ran2(&iseed);
                x3p = x3min + Lz*ran2(&iseed);
              }
              else
                x3p = x3l+pGrid->dx3/Npar*(kp+0.5);

              /*ath_pout(0, "[streaming3d]: x3ps: %f, %f, %f, %f.\n",
                                             x3l,
                                             pGrid->dx3,
                                             Npar,
                                             kp);*/

              pGrid->particle[p].property = 0;
              pGrid->particle[p].x1 = x1p;
              pGrid->particle[p].x2 = x2p;
              pGrid->particle[p].x3 = x3p;

              if ((ipert == 2) || (ipert == 3)) {
                pGrid->particle[p].x1 += -paramp*(cos(kx*x1p)*cos(-ky*x2p) -
                                          sin(kx*x1p)*sin(-ky*x2p));

                w1 = etavk * pert_XY(Rewx,Imwx,pGrid->particle[p].x1,x2p,t);
                w3 = etavk * pert_XY(Rewz,Imwz,pGrid->particle[p].x1,x2p,t);
                w2 = etavk * pert_XY(Rewy,Imwy,pGrid->particle[p].x1,x2p,t);
              } else {
                w1 = w2 = w3 = 0.0;
              }

              /*ath_pout(0, "[streaming3d]: particle pos: %f, %f, %f.\n",
                                             pGrid->particle[p].x1,
                                             pGrid->particle[p].x2,
                                             pGrid->particle[p].x3);*/

              pGrid->particle[p].v1 = wxNSH+w1;
              pGrid->particle[p].v2 = wyNSH+w2;

              pGrid->particle[p].v3 = w3;
#ifndef FARGO
              pGrid->particle[p].v2 -= qshear*Omega_0*x1p;
#endif
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


/* enroll gravitational potential function, shearing sheet BC functions */
  ShearingBoxPot = UnstratifiedDisk;
  ShBoxCoord = xy;

  if (myID_Comm_world == 0) {
    /* flush output file */
    sprintf(name, "%s_%d_%d.dat", "SelfgXY",Nx,ipert);
    FILE *fid = fopen(name,"w");
    fclose(fid);
#ifdef MPI_PARALLEL
    sprintf(name, "../%s_%d_%d.dat", "SelfgXY",Nx,ipert);
#endif
  }

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
  fwrite(&etavk, sizeof(Real),1,fp);
  fwrite(&uxNSH, sizeof(Real),1,fp); fwrite(&uyNSH, sizeof(Real),1,fp);
  fwrite(&wxNSH, sizeof(Real),1,fp); fwrite(&wyNSH, sizeof(Real),1,fp);
  if ((ipert==1) || (ipert==2)) {
    fwrite(&Reux, sizeof(Real),1,fp); fwrite(&Imux, sizeof(Real),1,fp);
    fwrite(&Reuy, sizeof(Real),1,fp); fwrite(&Imuy, sizeof(Real),1,fp);
    fwrite(&Reuz, sizeof(Real),1,fp); fwrite(&Imuz, sizeof(Real),1,fp);
    fwrite(&Rewx, sizeof(Real),1,fp); fwrite(&Imwx, sizeof(Real),1,fp);
    fwrite(&Rewy, sizeof(Real),1,fp); fwrite(&Imwy, sizeof(Real),1,fp);
    fwrite(&Rewz, sizeof(Real),1,fp); fwrite(&Imwz, sizeof(Real),1,fp);
    fwrite(&Rerho, sizeof(Real),1,fp);fwrite(&Imrho, sizeof(Real),1,fp);
    fwrite(&omg, sizeof(Real),1,fp);  fwrite(&s, sizeof(Real),1,fp);
    fwrite(&Iso_csound, sizeof(Real),1,fp);
    fwrite(&kx, sizeof(Real),1,fp);   fwrite(&ky, sizeof(Real),1,fp);
  }
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  ShearingBoxPot = UnstratifiedDisk;

  Omega_0 = par_getd("problem","omega");
  qshear = par_getd_def("problem","qshear",1.5);
  amp = par_getd_def("problem","amp",0.0);
  ipert = par_geti_def("problem","ipert", 1);

  Nx = pM->Nx[0]; /* used for particle selection function */
  Npar  = (int)(sqrt(par_geti("particle","parnumcell")));
  Npar3 = SQR(Npar)*Npar;
  downsamp = par_geti_def("problem","downsamp",Npar3);

  fread(&rho0, sizeof(Real),1,fp);  fread(&mratio, sizeof(Real),1,fp);
  fread(&etavk, sizeof(Real),1,fp);
  fread(&uxNSH, sizeof(Real),1,fp); fread(&uyNSH, sizeof(Real),1,fp);
  fread(&wxNSH, sizeof(Real),1,fp); fread(&wyNSH, sizeof(Real),1,fp);
  if ((ipert==1) || (ipert==2)) {
    fread(&Reux, sizeof(Real),1,fp); fread(&Imux, sizeof(Real),1,fp);
    fread(&Reuy, sizeof(Real),1,fp); fread(&Imuy, sizeof(Real),1,fp);
    fread(&Reuz, sizeof(Real),1,fp); fread(&Imuz, sizeof(Real),1,fp);
    fread(&Rewx, sizeof(Real),1,fp); fread(&Imwx, sizeof(Real),1,fp);
    fread(&Rewy, sizeof(Real),1,fp); fread(&Imwy, sizeof(Real),1,fp);
    fread(&Rewz, sizeof(Real),1,fp); fread(&Imwz, sizeof(Real),1,fp);
    fread(&Rerho, sizeof(Real),1,fp);fread(&Imrho, sizeof(Real),1,fp);
    fread(&omg, sizeof(Real),1,fp);  fread(&s, sizeof(Real),1,fp);
    fread(&Iso_csound, sizeof(Real),1,fp);  Iso_csound2 = SQR(Iso_csound);
    fread(&kx, sizeof(Real),1,fp);   fread(&ky, sizeof(Real),1,fp);
  }

  if (myID_Comm_world == 0)
#ifdef MPI_PARALLEL
    sprintf(name, "../%s_%d_%d.dat", "SelfgXY",Nx,ipert);
#else
    sprintf(name, "%s_%d_%d.dat", "SelfgXY",Nx,ipert);
#endif

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

/*! \fn static Real expr_rhopar(const GridS *pG, const int i, const int j, 
 *				   const int k)
 *  \brief dpar */
static Real expr_rhopar(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
return pG->Coup[k][j][i].grid_d;
}

/*! \fn static Real expr_dVx(const GridS *pG, const int i, const int j, 
 *			     const int k)
 *  \brief dVx */
static Real expr_dVx(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->U[k][j][i].M1/pG->U[k][j][i].d - uxNSH;
}

/*! \fn static Real expr_dVy(const GridS *pG, const int i, const int j, 
 *			     const int k)
 *  \brief dVy */
static Real expr_dVy(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M2/pG->U[k][j][i].d - uyNSH;
#else
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d -uyNSH + qshear*Omega_0*x1);
#endif
}

/*! \fn static Real expr_rhopardif(const GridS *pG, const int i, const int j, 
 *				   const int k)
 *  \brief difdpar */
static Real expr_rhopardif(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
return pG->Coup[k][j][i].grid_d - rho0*mratio;
}

/*! \fn static Real expr_dVxpar(const GridS *pG, const int i, const int j, 
 *				const int k)
 *  \brief dVxpar */
static Real expr_dVxpar(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return expr_V1par(pG,i,j,k) - wxNSH;
}

/*! \fn static Real expr_dVypar(const GridS *pG, const int i, const int j, 
 *			        const int k)
 *  \brief dVypar */
static Real expr_dVypar(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return expr_V2par(pG,i,j,k) - wyNSH;
#else
  return expr_V2par(pG,i,j,k) - wyNSH + qshear*Omega_0*x1;
#endif
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
  if(strcmp(expr,"difd")==0) return expr_rhodif;
  if(strcmp(expr,"dVx")==0) return expr_dVx;
  if(strcmp(expr,"dVy")==0) return expr_dVy;
  if(strcmp(expr,"difdpar")==0) return expr_rhopardif;
  if(strcmp(expr,"dVxpar")==0) return expr_dVxpar;
  if(strcmp(expr,"dVypar")==0) return expr_dVypar;
  if(strcmp(expr,"gphi1")==0) return expr_GradPhiX1;
  if(strcmp(expr,"gphi2")==0) return expr_GradPhiX2;
  if(strcmp(expr,"gphi3")==0) return expr_GradPhiX3;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  if(strcmp(name,"amp")==0) return OutputModeAmplitude;
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

/*--------------------------------------------------------------------------- */
/*! \fn void OutputModeAmplitude(Grid *pGrid, Domain *pDomain, OutputS *pOut)
 *  \brief output the perturbation amplitude */
void OutputModeAmplitude(MeshS *pM, OutputS *pOut)
{
  DomainS *pDomain = (DomainS*)&(pM->Domain[0][0]);
  GridS *pGrid = pDomain->Grid;

  FILE *fid;
  int i,j,k;
  Real dm,dparm,uxm,uym,uzm,wxm,wym,wzm;

  particle_to_grid(pDomain, property_all);

  dm=0.0; dparm=0.0; uxm=0.0; uym=0.0; uzm=0.0; wxm=0.0; wym=0.0; wzm=0.0;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
  for (i=pGrid->is; i<=pGrid->ie; i++) {
    dm  = MAX(dm,fabs(expr_rhodif(pGrid,i,j,k)));
    dparm=MAX(dparm,fabs(expr_rhopardif(pGrid,i,j,k)));
    uxm = MAX(uxm,fabs(expr_dVx(pGrid,i,j,k)));
    wxm = MAX(wxm,fabs(expr_dVxpar(pGrid,i,j,k)));
    uym = MAX(uym,fabs(expr_dVy(pGrid,i,j,k)));
    wym = MAX(wym,fabs(expr_dVypar(pGrid,i,j,k)));
    uzm = MAX(uzm,fabs(expr_V3(pGrid,i,j,k)));
    wzm = MAX(wzm,fabs(expr_V3par(pGrid,i,j,k)));
  }}}

#ifdef MPI_PARALLEL
  Real sendbuf[8], recvbuf[8];
  int err;
  sendbuf[0]=dm;	sendbuf[1]=dparm;
  sendbuf[2]=uxm;	sendbuf[3]=wxm;
  sendbuf[4]=uym;	sendbuf[5]=wym;
  sendbuf[6]=uzm;	sendbuf[7]=wzm;

  err = MPI_Reduce(sendbuf,recvbuf,8,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(err) ath_error("[streaming3d]: MPI_Reduce returned error code %d\n",err);

  if (myID_Comm_world == 0) {
    dm=recvbuf[0];	dparm=recvbuf[1];
    uxm=recvbuf[2];	wxm=recvbuf[3];
    uym=recvbuf[4];	wym=recvbuf[5];
    uzm=recvbuf[6];	wzm=recvbuf[7];
  }
#endif

  if (myID_Comm_world  == 0) {
    fid = fopen(name,"a+");
    fprintf(fid,"%e %e %e %e %e %e %e %e %e\n",
                pGrid->time,dm,dparm,uxm,wxm,uym,wym,uzm,wzm);
    fclose(fid);
  }

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
