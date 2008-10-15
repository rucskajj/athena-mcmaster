#include "../copyright.h"
/*==============================================================================
 * FILE: resistivity.c
 *
 * PURPOSE: Implements explicit Ohmic resistivity, that is
 *      dB/dt = -Curl(\eta J)    where J=Curl(B).
 *   Functions are called by integrate_diffusion() in the main loop, which
 *   coordinates adding all diffusion operators (viscosity, resistivity, thermal
 *   conduction) using operator splitting.
 *
 *   An explicit timestep limit must be applied if these routines are used.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  ohmic_resistivity_1d()
 *  ohmic_resistivity_2d()
 *  ohmic_resistivity_3d()
 *  ohmic_resistivity_init() - allocates memory needed
 *  ohmic_resistivity_destruct() - frees memory used
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef OHMIC
#ifdef ADIABATIC
#error : resistivity only works for isothermal EOS.
#endif /* ADIABATIC */
#ifdef HYDRO
#error : Ohmic resistivity only works for MHD.
#endif /* HYDRO */
#endif

/* The resistive emfs, contained in special structure */
typedef struct ThreeDVect_t{
  Real x;
  Real y;
  Real z;
}ThreeDVect;
static ThreeDVect ***emf=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_1d:
 */

void ohmic_resistivity_1d(Grid *pG, Domain *pD)
{
#ifdef OHMIC
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1;

/*--- Step 1 -------------------------------------------------------------------
 * Compute resistive EMF.  Note:
 *   emf.x = eta*J1 = 0
 *   emf.y = eta*J2 = eta_R*(-dB3/dx1)
 *   emf.z = eta*J3 = eta_R*(dB2/dx1)
 * Note emf.y and emf.z use B3c and B2c respectively
 */

  for (i=is; i<=ie+1; i++) {
    emf[ks][js][i].x = 0.0;
    emf[ks][js][i].y = -(pG->U[ks][js][i].B3c - pG->U[ks][js][i-1].B3c)/pG->dx1;
    emf[ks][js][i].z =  (pG->U[ks][js][i].B2c - pG->U[ks][js][i-1].B2c)/pG->dx1;
/* Multiple components by constant \eta_R */
    emf[ks][js][i].y *= eta_R;
    emf[ks][js][i].z *= eta_R;
  }

/*--- Step 2 -------------------------------------------------------------------
 * CT update of magnetic field using resistive EMF.  In 1D, this reduces to
 * centered differences for the resistive fluxes of B2c and B3c
 */

  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].B2c += dtodx1*(emf[ks][js][i+1].z - emf[ks][js][i].z);
    pG->U[ks][js][i].B3c -= dtodx1*(emf[ks][js][i+1].y - emf[ks][js][i].y);
  }

/*--- Step 3 -------------------------------------------------------------------
 * For consistency, set B2i and B3i to cell-centered values. */

  for (i=is; i<=ie; i++) {
    pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
    pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;
  }
#endif /* OHMIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_2d:
 */

void ohmic_resistivity_2d(Grid *pG, Domain *pD)
{
#ifdef OHMIC
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;

/*--- Step 1 -------------------------------------------------------------------
 * Compute resistive EMF.  Note:
 *   emf.x = eta*J1 = eta_R*(dB3/dx2)
 *   emf.y = eta*J2 = eta_R*(-dB3/dx1)
 *   emf.z = eta*J3 = eta_R*(dB2/dx1 - dB1/dx2)
 * Note emf.x and emf.y use B3c
 */

  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      emf[ks][j][i].x =  (pG->U[ks][j][i].B3c - pG->U[ks][j-1][i].B3c)/pG->dx2;
      emf[ks][j][i].y = -(pG->U[ks][j][i].B3c - pG->U[ks][j][i-1].B3c)/pG->dx1;

      emf[ks][j][i].z = (pG->B2i[ks][j][i] - pG->B2i[ks][j  ][i-1])/pG->dx1 -
                        (pG->B1i[ks][j][i] - pG->B1i[ks][j-1][i  ])/pG->dx2;
/* Multiple components by constant \eta_R */
      emf[ks][j][i].x *= eta_R;
      emf[ks][j][i].y *= eta_R;
      emf[ks][j][i].z *= eta_R;
    }
  }

/*--- Step 2 -------------------------------------------------------------------
 * CT update of magnetic field using resistive EMF.  This is identical to the
 * CT update in the 2D integrators: dB/dt = -Curl(E).  For B3, the CT formula
 * reduces to centered differences for the diffusive (resistive) fluxes of B3c.
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B1i[ks][j][i] -= dtodx2*(emf[ks][j+1][i  ].z - emf[ks][j][i].z);
      pG->B2i[ks][j][i] += dtodx1*(emf[ks][j  ][i+1].z - emf[ks][j][i].z);

      pG->U[ks][j][i].B3c += dtodx2*(emf[ks][j+1][i  ].x - emf[ks][j][i].x) -
                             dtodx1*(emf[ks][j  ][i+1].y - emf[ks][j][i].y);
    }
    pG->B1i[ks][j][ie+1] -= dtodx2*(emf[ks][j+1][ie+1].z - emf[ks][j][ie+1].z);
  }
  for (i=is; i<=ie; i++) {
    pG->B2i[ks][je+1][i] += dtodx1*(emf[ks][je+1][i+1].z - emf[ks][je+1][i].z);
  }

/*--- Step 3 -------------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].B1c = 0.5*(pG->B1i[ks][j][i] + pG->B1i[ks][j][i+1]);
      pG->U[ks][j][i].B2c = 0.5*(pG->B2i[ks][j][i] + pG->B2i[ks][j+1][i]);
/* Set the 3-interface magnetic field equal to the cell center field. */
      pG->B3i[ks][j][i] = pG->U[ks][j][i].B3c;
    }
  }
#endif /* OHMIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_3d:
 */

void ohmic_resistivity_3d(Grid *pG, Domain *pD)
{
#ifdef OHMIC
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;

/*--- Step 1 -------------------------------------------------------------------
 * Compute resistive EMFs.  Note:
 *   emf.x = eta*J1 = eta_R*(dB3/dx2 - dB2/dx3)
 *   emf.y = eta*J2 = eta_R*(dB1/dx3 - dB3/dx1)
 *   emf.z = eta*J3 = eta_R*(dB2/dx1 - dB1/dx2)
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        emf[k][j][i].x = (pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ])/pG->dx2 -
                         (pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ])/pG->dx3;
        emf[k][j][i].y = (pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ])/pG->dx3 -
                         (pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1])/pG->dx1;
        emf[k][j][i].z = (pG->B2i[k][j][i] - pG->B2i[k  ][j  ][i-1])/pG->dx1 -
                         (pG->B1i[k][j][i] - pG->B1i[k  ][j-1][i  ])/pG->dx2;
/* Multiple components by constant \eta_R */
        emf[k][j][i].x *= eta_R;
        emf[k][j][i].y *= eta_R;
        emf[k][j][i].z *= eta_R;
      }
    }
  }

/*--- Step 2 -------------------------------------------------------------------
 * CT update of magnetic field using resistive EMFs.  This is identical to the
 * CT update in the integrators: dB/dt = -Curl(E)
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->B1i[k][j][i] += dtodx3*(emf[k+1][j  ][i  ].y - emf[k][j][i].y) -
                            dtodx2*(emf[k  ][j+1][i  ].z - emf[k][j][i].z);
        pG->B2i[k][j][i] += dtodx1*(emf[k  ][j  ][i+1].z - emf[k][j][i].z) -
                            dtodx3*(emf[k+1][j  ][i  ].x - emf[k][j][i].x);
        pG->B3i[k][j][i] += dtodx2*(emf[k  ][j+1][i  ].x - emf[k][j][i].x) -
                            dtodx1*(emf[k  ][j  ][i+1].y - emf[k][j][i].y);
      }
      pG->B1i[k][j][ie+1] +=
        dtodx3*(emf[k+1][j  ][ie+1].y - emf[k][j][ie+1].y) -
        dtodx2*(emf[k  ][j+1][ie+1].z - emf[k][j][ie+1].z);
    }
    for (i=is; i<=ie; i++) {
      pG->B2i[k][je+1][i] +=
        dtodx1*(emf[k  ][je+1][i+1].z - emf[k][je+1][i].z) -
        dtodx3*(emf[k+1][je+1][i  ].x - emf[k][je+1][i].x);
    }
  }
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B3i[ke+1][j][i] +=
        dtodx2*(emf[ke+1][j+1][i  ].x - emf[ke+1][j][i].x) -
        dtodx1*(emf[ke+1][j  ][i+1].y - emf[ke+1][j][i].y);
    }
  }

/*--- Step 3 -------------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
      }
    }
  }
#endif /* OHMIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_init: Allocate temporary arrays
 */

void ohmic_resistivity_init(int nx1, int nx2, int nx3)
{
#ifdef OHMIC
  int Nx1 = nx1 + 2;
  int Nx2, Nx3;

  if (nx2 > 1) Nx2 = nx2 + 2;
  if (nx3 > 1) Nx3 = nx3 + 2;
  
  if ((emf = (ThreeDVect***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ThreeDVect)))
    == NULL) goto on_error;
  return;

  on_error:
  ohmic_resistivity_destruct();
  ath_error("[ohmic_resisticvity_init]: malloc returned a NULL pointer\n");
#endif /* OHMIC */
  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_destruct: Free temporary arrays
 */

void ohmic_resistivity_destruct(void)
{
#ifdef OHMIC
  if (emf != NULL) free_3d_array(emf);
#endif /* OHMIC */
  return;
}