#include "../copyright.h"
/*============================================================================*/
/*! \file bvals_grav.c
 *  \brief Sets boundary conditions (quantities in ghost zones) for the
 *   gravitational potential on each edge of a Grid. 
 *
 * PURPOSE: Sets boundary conditions (quantities in ghost zones) for the
 *   gravitational potential on each edge of a Grid.  See comments at
 *   start of bvals_mhd.c for more details.
 * The only BC functions implemented here are for:
 *- 1 = reflecting, 4 = periodic, and MPI boundaries
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - bvals_grav()      - calls appropriate functions to set ghost cells
 * - bvals_grav_init() - sets function pointers used by bvals_grav()
 * - bvals_grav_fun()  - enrolls a pointer to a user-defined BC function */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "prototypes.h"
#include "../prototypes.h"
#include "../globals.h"

/* The functions in this file will only work with SELF_GRAVITY */
#ifdef SELF_GRAVITY

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
#endif /* MPI_PARALLEL */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   calcGradPhi() - calculate GradPhi for particle gravity
 *   reflect_Phi_???()  - apply reflecting BCs at boundary ???
 *   periodic_Phi_???() - apply periodic BCs at boundary ???
 *   pack_Phi_???()   - pack data for MPI non-blocking send at ??? boundary
 *   unpack_Phi_???() - unpack data for MPI non-blocking receive at ??? boundary
 *============================================================================*/
#ifdef PARTICLES
static void calcGradPhi(GridS *pGrid);
static void copyPhiPar(GridS *pGrid);
#endif

static void reflect_Phi_ix1(GridS *pG);
static void reflect_Phi_ox1(GridS *pG);
static void reflect_Phi_ix2(GridS *pG);
static void reflect_Phi_ox2(GridS *pG);
static void reflect_Phi_ix3(GridS *pG);
static void reflect_Phi_ox3(GridS *pG);

static void periodic_Phi_ix1(GridS *pG);
static void periodic_Phi_ox1(GridS *pG);
static void periodic_Phi_ix2(GridS *pG);
static void periodic_Phi_ox2(GridS *pG);
static void periodic_Phi_ix3(GridS *pG);
static void periodic_Phi_ox3(GridS *pG);

static void obc_fft_Phi_ix1(GridS *pG);
static void obc_fft_Phi_ox1(GridS *pG);
static void obc_fft_Phi_ix2(GridS *pG);
static void obc_fft_Phi_ox2(GridS *pG);
static void obc_fft_Phi_ix3(GridS *pG);
static void obc_fft_Phi_ox3(GridS *pG);

static void ProlongateLater(GridS *pG);

#ifdef MPI_PARALLEL
static void pack_Phi_ix1(GridS *pG);
static void pack_Phi_ox1(GridS *pG);
static void pack_Phi_ix2(GridS *pG);
static void pack_Phi_ox2(GridS *pG);
static void pack_Phi_ix3(GridS *pG);
static void pack_Phi_ox3(GridS *pG);

static void unpack_Phi_ix1(GridS *pG);
static void unpack_Phi_ox1(GridS *pG);
static void unpack_Phi_ix2(GridS *pG);
static void unpack_Phi_ox2(GridS *pG);
static void unpack_Phi_ix3(GridS *pG);
static void unpack_Phi_ox3(GridS *pG);
#endif /* MPI_PARALLEL */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void bvals_grav(DomainS *pD)
 *  \brief Calls appropriate functions to set ghost zones.  
 *
 *   The function
 *   pointers (*_GBCFun) are set during initialization by bvals_grav_init() to
 *   be one of the functions corresponding to reflecting or periodic.  If the
 *   left- or right-Grid ID numbers are >= 1 (neighboring grids exist), then MPI
 *   calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void bvals_grav(DomainS *pD)
{
  GridS *pGrid = (pD->Grid);
#ifdef SHEARING_BOX
  int myL,myM,myN;
#endif
#ifdef MPI_PARALLEL
  int cnt1, cnt2, cnt3, cnt, ierr, mIndex;
#endif /* MPI_PARALLEL */

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pGrid->Nx[0] > 1){

#ifdef MPI_PARALLEL

    cnt = nghost*(pGrid->Nx[1])*(pGrid->Nx[2]);

/* MPI blocks to both left and right */
    if (pGrid->rx1_Gid >= 0 && pGrid->lx1_Gid >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_Phi_ix1(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_Phi_ox1(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix1(pGrid);
      if (mIndex == 1) unpack_Phi_ox1(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix1(pGrid);
      if (mIndex == 1) unpack_Phi_ox1(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx1_Gid >= 0 && pGrid->lx1_Gid < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_Phi_ox1(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix1_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_Phi_ox1(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx1_Gid < 0 && pGrid->lx1_Gid >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_Phi_ix1(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox1_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_Phi_ix1(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx1_Gid < 0 && pGrid->lx1_Gid < 0) {
      (*(pD->ix1_GBCFun))(pGrid);
      (*(pD->ox1_GBCFun))(pGrid);
    }

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pGrid->Nx[1] > 1){

#ifdef MPI_PARALLEL

    cnt = (pGrid->Nx[0] + 2*nghost)*nghost*(pGrid->Nx[2]);

/* MPI blocks to both left and right */
    if (pGrid->rx2_Gid >= 0 && pGrid->lx2_Gid >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_Phi_ix2(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_Phi_ox2(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix2(pGrid);
      if (mIndex == 1) unpack_Phi_ox2(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix2(pGrid);
      if (mIndex == 1) unpack_Phi_ox2(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx2_Gid >= 0 && pGrid->lx2_Gid < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_Phi_ox2(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix2_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_Phi_ox2(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx2_Gid < 0 && pGrid->lx2_Gid >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_Phi_ix2(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox2_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_Phi_ix2(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx2_Gid < 0 && pGrid->lx2_Gid < 0) {
      (*(pD->ix2_GBCFun))(pGrid);
      (*(pD->ox2_GBCFun))(pGrid);
    }

/* shearing sheet BCs; function defined in problem generator */
#ifdef SHEARING_BOX
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
    if (myL == 0) {
      ShearingSheet_grav_ix1(pD);
    }
    if (myL == (pD->NGrid[0]-1)) {
      ShearingSheet_grav_ox1(pD);
    }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pGrid->Nx[2] > 1){

#ifdef MPI_PARALLEL

    cnt = (pGrid->Nx[0] + 2*nghost)*(pGrid->Nx[1] + 2*nghost)*nghost;

/* MPI blocks to both left and right */
    if (pGrid->rx3_Gid >= 0 && pGrid->lx3_Gid >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_Phi_ix3(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_Phi_ox3(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix3(pGrid);
      if (mIndex == 1) unpack_Phi_ox3(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_Phi_ix3(pGrid);
      if (mIndex == 1) unpack_Phi_ox3(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx3_Gid >= 0 && pGrid->lx3_Gid < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_Gid,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_Phi_ox3(pGrid);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_Gid,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix3_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_Phi_ox3(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx3_Gid < 0 && pGrid->lx3_Gid >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_Gid,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_Phi_ix3(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_Gid,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox3_GBCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_Phi_ix3(pGrid);
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx3_Gid < 0 && pGrid->lx3_Gid < 0) {
      (*(pD->ix3_GBCFun))(pGrid);
      (*(pD->ox3_GBCFun))(pGrid);
    }

  }

  
#ifdef PARTICLES
  calcGradPhi(pGrid);
  copyPhiPar(pGrid);
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_grav_init(MeshS *pM) 
 *  \brief Sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI
 */

void bvals_grav_init(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int i,nl,nd,irefine;
#ifdef MPI_PARALLEL
  int myL,myM,myN,l,m,n,nx1t,nx2t,nx3t,size;
  int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
#endif /* MPI_PARALLEL */

/* Cycle through all the Domains that have active Grids on this proc */

  for (nl=0; nl<(pM->NLevels); nl++){
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
  if (pM->Domain[nl][nd].Grid != NULL) {
    pD = (DomainS*)&(pM->Domain[nl][nd]);  /* ptr to Domain */
    pG = pM->Domain[nl][nd].Grid;          /* ptr to Grid */
    irefine = 1;
    for (i=1;i<=nl;i++) irefine *= 2;   /* C pow fn only takes doubles !! */
#ifdef MPI_PARALLEL
/* get (l,m,n) coordinates of Grid being updated on this processor */
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
#endif /* MPI_PARALLEL */

/* Set function pointers for physical boundaries in x1-direction */

    if(pG->Nx[0] > 1) {

/*---- ix1 boundary ----------------------------------------------------------*/

      if(pD->ix1_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if(pD->Disp[0] != 0) {
          pD->ix1_GBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain */
        } else {
#ifndef SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ix1){

          case 1: /* Reflecting */
            pD->ix1_GBCFun = reflect_Phi_ix1;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ix1 = 2 not implemented\n");
          break;

          case 4: /* Periodic */
            pD->ix1_GBCFun = periodic_Phi_ix1;
#ifdef MPI_PARALLEL
            if(pG->lx1_Gid < 0 && pD->NGrid[0] > 1){
              pG->lx1_Gid = pD->GData[myN][myM][pD->NGrid[0]-1].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ix1_GBCFun = reflect_Phi_ix1;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ix1 = %d unknown\n",
            pM->BCFlag_ix1);
          }
#else /* SELF_GRAVITY_USING_FFT_OBC */
        pD->ix1_GBCFun = obc_fft_Phi_ix1;
#endif
        }
      }

/*---- ox1 boundary ----------------------------------------------------------*/

      if(pD->ox1_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) {
          pD->ox1_GBCFun = ProlongateLater;

/* Domain is at R-edge of root Domain */
        } else {
#ifndef SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ox1){

          case 1: /* Reflecting */
            pD->ox1_GBCFun = reflect_Phi_ox1;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ox1 = 2 not implemented\n");
          break;

          case 4: /* Periodic */
            pD->ox1_GBCFun = periodic_Phi_ox1;
#ifdef MPI_PARALLEL
            if(pG->rx1_Gid < 0 && pD->NGrid[0] > 1){
              pG->rx1_Gid = pD->GData[myN][myM][0].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ox1_GBCFun = reflect_Phi_ox1;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ox1 = %d unknown\n",
            pM->BCFlag_ox1);
          }
#else /* SELF_GRAVITY_USING_FFT_OBC */
          pD->ox1_GBCFun = obc_fft_Phi_ox1;
#endif
        }
      }
    }

/* Set function pointers for physical boundaries in x2-direction */

    if(pG->Nx[1] > 1) {

/*---- ix2 boundary ----------------------------------------------------------*/

      if(pD->ix2_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if(pD->Disp[1] != 0) {
          pD->ix2_GBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain */
        } else {
#ifndef SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ix2){

          case 1: /* Reflecting */
            pD->ix2_GBCFun = reflect_Phi_ix2;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ix2 = 2 not implemented\n");
          break;

          case 4: /* Periodic */
            pD->ix2_GBCFun = periodic_Phi_ix2;
#ifdef MPI_PARALLEL
            if(pG->lx2_Gid < 0 && pD->NGrid[1] > 1){
              pG->lx2_Gid = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ix2_GBCFun = reflect_Phi_ix2;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ix2 = %d unknown\n",
            pM->BCFlag_ix2);
          }
#else /* SELF_GRAVITY_USING_FFT_OBC */
          pD->ix2_GBCFun = obc_fft_Phi_ix2;
#endif
        }
      }

/*---- ox2 boundary ----------------------------------------------------------*/

    if(pD->ox2_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) {
          pD->ox2_GBCFun = ProlongateLater;

/* Domain is at R-edge of root Domain */
        } else {
#ifndef SELF_GRAVITY_USING_FFT_OBC
          switch(pM->BCFlag_ox2){

          case 1: /* Reflecting */
            pD->ox2_GBCFun = reflect_Phi_ox2;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ox2 = 2 not implemented\n");
          break;

          case 4: /* Periodic */
            pD->ox2_GBCFun = periodic_Phi_ox2;
#ifdef MPI_PARALLEL
            if(pG->rx2_Gid < 0 && pD->NGrid[1] > 1){
              pG->rx2_Gid = pD->GData[myN][0][myL].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ox2_GBCFun = reflect_Phi_ox2;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ox2 = %d unknown\n",
            pM->BCFlag_ox2);
          }
#else /* SELF_GRAVITY_USING_FFT_OBC */
          pD->ox2_GBCFun = obc_fft_Phi_ox2;
#endif
        }
      }
    }

/* Set function pointers for physical boundaries in x3-direction */

    if(pG->Nx[2] > 1) {

/*---- ix3 boundary ----------------------------------------------------------*/

      if(pD->ix3_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if(pD->Disp[2] != 0) {
          pD->ix3_GBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain */
        } else {
#if defined(SELF_GRAVITY_USING_FFT_OBC) || defined(SELF_GRAVITY_USING_FFT_DISK)
          pD->ix3_GBCFun = obc_fft_Phi_ix3;
#else
          switch(pM->BCFlag_ix3){

          case 1: /* Reflecting */
            pD->ix3_GBCFun = reflect_Phi_ix3;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ix3 = 2 not defined\n");
          break;

          case 4: /* Periodic */
            pD->ix3_GBCFun = periodic_Phi_ix3;
#ifdef MPI_PARALLEL
            if(pG->lx3_Gid < 0 && pD->NGrid[2] > 1){
              pG->lx3_Gid = pD->GData[pD->NGrid[2]-1][myM][myL].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ix3_GBCFun = reflect_Phi_ix3;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ix3 = %d unknown\n",
            pM->BCFlag_ix3);
          }
#endif /* SELF_GRAVITY_USING_FFT_OBC || SELF_GRAVITY_USING_FFT_DISK*/
        }
      }

/*---- ox3 boundary ----------------------------------------------------------*/

    if(pD->ox3_GBCFun == NULL){

/* Domain boundary is in interior of root */
        if((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) {
          pD->ox3_GBCFun = ProlongateLater;

/* Domain is at R-edge of root Domain */
        } else {
#if defined(SELF_GRAVITY_USING_FFT_OBC) || defined(SELF_GRAVITY_USING_FFT_DISK)
          pD->ox3_GBCFun = obc_fft_Phi_ox3;
#else
          switch(pM->BCFlag_ox3){

          case 1: /* Reflecting */
            pD->ox3_GBCFun = reflect_Phi_ox3;
          break;

          case 2: /* Outflow */
            ath_error("[bvals_grav_init]: BCFlag_ox3 = 2 not defined\n");
          break;

          case 4: /* Periodic */
            pD->ox3_GBCFun = periodic_Phi_ox3;
#ifdef MPI_PARALLEL
            if(pG->rx3_Gid < 0 && pD->NGrid[2] > 1){
              pG->rx3_Gid = pD->GData[0][myM][myL].ID_Comm_Domain;
            }
#endif /* MPI_PARALLEL */
          break;

          case 5: /* Reflecting, B_normal!=0 */
            pD->ox3_GBCFun = reflect_Phi_ox3;
          break;

          default:
            ath_error("[bvals_grav_init]: BCFlag_ox3 = %d unknown\n",
            pM->BCFlag_ox3);
          }
#endif /* SELF_GRAVITY_USING_FFT_OBC || SELF_GRAVITY_USING_FFT_DISK*/
        }
      }
    }

/* Figure out largest size needed for send/receive buffers with MPI ----------*/

#ifdef MPI_PARALLEL

    for (n=0; n<(pD->NGrid[2]); n++){
    for (m=0; m<(pD->NGrid[1]); m++){
      for (l=0; l<(pD->NGrid[0]); l++){

/* x1cnt is surface area of x1 faces */
        if(pD->NGrid[0] > 1){
          nx2t = pD->GData[n][m][l].Nx[1];
          if(nx2t > 1) nx2t += 1;

          nx3t = pD->GData[n][m][l].Nx[2];
          if(nx3t > 1) nx3t += 1;

          if(nx2t*nx3t > x1cnt) x1cnt = nx2t*nx3t;
        }

/* x2cnt is surface area of x2 faces */
        if(pD->NGrid[1] > 1){
          nx1t = pD->GData[n][m][l].Nx[0];
          if(nx1t > 1) nx1t += 2*nghost;

          nx3t = pD->GData[n][m][l].Nx[2];
          if(nx3t > 1) nx3t += 1;

          if(nx1t*nx3t > x2cnt) x2cnt = nx1t*nx3t;
        }

/* x3cnt is surface area of x3 faces */
        if(pD->NGrid[2] > 1){
          nx1t = pD->GData[n][m][l].Nx[0];
          if(nx1t > 1) nx1t += 2*nghost;

          nx2t = pD->GData[n][m][l].Nx[1];
          if(nx2t > 1) nx2t += 2*nghost;

          if(nx1t*nx2t > x3cnt) x3cnt = nx1t*nx2t;
        }
      }
    }}
#endif /* MPI_PARALLEL */

  }}}  /* End loop over all Domains with active Grids -----------------------*/

#ifdef MPI_PARALLEL
/* Allocate memory for send/receive buffers and MPI_Requests */

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  size *= nghost; /* Multiply by the third dimension */

  if (size > 0) {
    if((send_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate recv buffer\n");
  }

  if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate send MPI_Request array\n");

#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_grav_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
 *  \brief Sets function pointers for user-defined BCs in problem 
 */

void bvals_grav_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
{
  switch(dir){
  case left_x1:
    pD->ix1_GBCFun = prob_bc;
    break;
  case right_x1:
    pD->ox1_GBCFun = prob_bc;
    break;
  case left_x2:
    pD->ix2_GBCFun = prob_bc;
    break;
  case right_x2:
    pD->ox2_GBCFun = prob_bc;
    break;
  case left_x3:
    pD->ix3_GBCFun = prob_bc;
    break;
  case right_x3:
    pD->ox3_GBCFun = prob_bc;
    break;
  default:
    ath_perr(-1,"[bvals_grav_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   calcGradPhi() - calculate GradPhi for particle gravity
 *   reflecting_???
 *   periodic_???
 *   send_???
 *   receive_???
 * where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 */

#ifdef PARTICLES
/*----------------------------------------------------------------------------*/
/*! \fn static void calcGradPhi(GridS *pGrid)
 *  \brief Calculate GradPhi if particle module is compiled
 */
static void calcGradPhi(GridS *pGrid)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;

  Real3Vect cell1;              /* one over dx1, dx2, dx3 */
  /* cell1 is a shortcut expressions as well as dimension indicator */
  if (pGrid->Nx[0] > 1)  cell1.x1 = 1.0/pGrid->dx1;  else cell1.x1 = 0.0;
  if (pGrid->Nx[1] > 1)  cell1.x2 = 1.0/pGrid->dx2;  else cell1.x2 = 0.0;
  if (pGrid->Nx[2] > 1)  cell1.x3 = 1.0/pGrid->dx3;  else cell1.x3 = 0.0;


  for (k=ks-(nghost-1); k<=ke+(nghost-1); k++){
    for (j=js-(nghost-1); j<=je+(nghost-1); j++){
      for (i=is-(nghost-1); i<=ie+(nghost-1); i++){

        //if(pGrid->time == 0.0){ /* First step, no Phi_old for average */
          pGrid->GradPhiX1[k][j][i] = -0.5*cell1.x1* \
           ( pGrid->Phi[k][j][i+1] - pGrid->Phi[k][j][i-1] );

          pGrid->GradPhiX2[k][j][i] = -0.5*cell1.x2* \
           ( pGrid->Phi[k][j+1][i] - pGrid->Phi[k][j-1][i] );

          pGrid->GradPhiX3[k][j][i] = -0.5*cell1.x3* \
           ( pGrid->Phi[k+1][j][i] - pGrid->Phi[k-1][j][i] );
        //}
        /*else{ Average Phi with Phi_old 
          pGrid->GradPhiX1[k][j][i] = -0.5*cell1.x1* \
            ( 0.5*(pGrid->Phi[k][j][i+1]+pGrid->Phi_old[k][j][i+1]) - \
              0.5*(pGrid->Phi[k][j][i-1]+pGrid->Phi_old[k][j][i-1]));

          pGrid->GradPhiX2[k][j][i] = -0.5*cell1.x2* \
            ( 0.5*(pGrid->Phi[k][j+1][i]+pGrid->Phi_old[k][j+1][i]) - \
              0.5*(pGrid->Phi[k][j-1][i]+pGrid->Phi_old[k][j-1][i]));

          pGrid->GradPhiX3[k][j][i] = -0.5*cell1.x3* \
            ( 0.5*(pGrid->Phi[k+1][j][i]+pGrid->Phi_old[k+1][j][i]) - \
              0.5*(pGrid->Phi[k-1][j][i]+pGrid->Phi_old[k-1][j][i]));

        if(k == 10){
          ath_pout(0,"[i,j,k]: [%d,%d,%d] ; Phi, Phi_old: %g %g ; cell1.x1: %g; dx1: %g\n", i,j,k, pGrid->Phi[k][j][i], pGrid->Phi_old[k][j][i], cell1.x1, pGrid->dx1);
        }


        }*/
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn static void calcPhiPar(GridS *pGrid)
 *  \brief Copy Phi to Phi_par to remove gravity force calculation for gas
 */
static void copyPhiPar(GridS *pGrid)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;

  for (k=ks-(nghost); k<=ke+(nghost); k++){
    for (j=js-(nghost); j<=je+(nghost); j++){
      for (i=is-(nghost); i<=ie+(nghost); i++){
        pGrid->Phi_par[k][j][i] = pGrid->Phi[k][j][i];

        pGrid->Phi[k][j][i] = 0.0;
        pGrid->Phi_old[k][j][i] = 0.0;
      }
    }
  }
}
#endif /* PARTICLES */

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ix1(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x1 boundary (ibc_x1=1)
 */

static void reflect_Phi_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->Phi[k][j][is-i] = pGrid->Phi[k][j][is+(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ox1(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x1 boundary (obc_x1=1)
 */

static void reflect_Phi_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->Phi[k][j][ie+i] = pGrid->Phi[k][j][ie-(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ix2(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x2 boundary (ibc_x2=1)
 */

static void reflect_Phi_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->Phi[k][js-j][i]    =  pGrid->Phi[k][js+(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ox2(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x2 boundary (obc_x2=1)
 */

static void reflect_Phi_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->Phi[k][je+j][i] = pGrid->Phi[k][je-(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ix3(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Inner x3 boundary (ibc_x3=1)
 */

static void reflect_Phi_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->Phi[ks-k][j][i] = pGrid->Phi[ks+(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_Phi_ox3(GridS *pGrid)
 *  \brief REFLECTING boundary conditions, Outer x3 boundary (obc_x3=1)
 */

static void reflect_Phi_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->Phi[ke+k][j][i] = pGrid->Phi[ke-(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ix1(GridS *pGrid)
 *  \brief PERIODIC boundary conditions, Inner x1 boundary (ibc_x1=4)
 */

static void periodic_Phi_ix1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->Phi[k][j][is-i] = pGrid->Phi[k][j][ie-(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ox1(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x1 boundary (obc_x1=4)
 */

static void periodic_Phi_ox1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->Phi[k][j][ie+i] = pGrid->Phi[k][j][is+(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ix2(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Inner x2 boundary (ibc_x2=4)
 */

static void periodic_Phi_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->Phi[k][js-j][i] = pGrid->Phi[k][je-(j-1)][i];
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ox2(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x2 boundary (obc_x2=4)
 */

static void periodic_Phi_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->Phi[k][je+j][i] = pGrid->Phi[k][js+(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ix3(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Inner x3 boundary (ibc_x3=4)
 */

static void periodic_Phi_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->Phi[ks-k][j][i] = pGrid->Phi[ke-(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_Phi_ox3(GridS *pGrid)
 *  \brief PERIODIC boundary conditions (cont), Outer x3 boundary (obc_x3=4)
 */

static void periodic_Phi_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->Phi[ke+k][j][i] = pGrid->Phi[ks+(k-1)][j][i];
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------
 * OPEN BOUNDARY CONDITION FUNCTIONS
 *----------------------------------------------------------------------------*/

/* For linear extrapolation, just comment out the following line */
#define QUADRATIC_EXTRAPOLATION

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Inner x1 boundary
 *
 * This BC uses a linear (quadratic) extrapolation of Phi normal to each face.
 * NOTE:  This version requires AT LEAST 2(3) active zones in each direction! */
static void obc_fft_Phi_ix1(GridS *pG)
{
  int i,is=pG->is;
  int j,js=pG->js,je=pG->je;
  int k,ks=pG->ks,ke=pG->ke;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->Phi[k][j][is-i] = 3.0*pG->Phi[k][j][is-i+1] - 3.0*pG->Phi[k][j][is-i+2] + pG->Phi[k][j][is-i+3];
#else  /* LINEAR EXTRAPOLATION */
        pG->Phi[k][j][is-i] = 2.0*pG->Phi[k][j][is-i+1] - pG->Phi[k][j][is-i+2];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Outer x1 boundary
 */
static void obc_fft_Phi_ox1(GridS *pG)
{
  int i,ie=pG->ie;
  int j,js=pG->js,je=pG->je;
  int k,ks=pG->ks,ke=pG->ke;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->Phi[k][j][ie+i] = 3.0*pG->Phi[k][j][ie+i-1] - 3.0*pG->Phi[k][j][ie+i-2] + pG->Phi[k][j][ie+i-3];
#else  /* LINEAR EXTRAPOLATION */
        pG->Phi[k][j][ie+i] = 2.0*pG->Phi[k][j][ie+i-1] - pG->Phi[k][j][ie+i-2];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Innter x2 boundary
 */
static void obc_fft_Phi_ix2(GridS *pG)
{
  int i,is=pG->is,ie=pG->ie;
  int j,js=pG->js;
  int k,ks=pG->ks,ke=pG->ke;
  
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->Phi[k][js-j][i] = 3.0*pG->Phi[k][js-j+1][i] - 3.0*pG->Phi[k][js-j+2][i] + pG->Phi[k][js-j+3][i];
#else  /* LINEAR EXTRAPOLATION */
        pG->Phi[k][js-j][i] = 2.0*pG->Phi[k][js-j+1][i] - pG->Phi[k][js-j+2][i];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Outer x2 boundary
 */
static void obc_fft_Phi_ox2(GridS *pG)
{
  int i,is=pG->is,ie=pG->ie;
  int j,je=pG->je;
  int k,ks=pG->ks,ke=pG->ke;
  
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->Phi[k][je+j][i] = 3.0*pG->Phi[k][je+j-1][i] - 3.0*pG->Phi[k][je+j-2][i] + pG->Phi[k][je+j-3][i];
#else  /* LINEAR EXTRAPOLATION */
        pG->Phi[k][je+j][i] = 2.0*pG->Phi[k][je+j-1][i] - pG->Phi[k][je+j-2][i];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Inner x3 boundary
 */
static void obc_fft_Phi_ix3(GridS *pG)
{
  int i,is=pG->is,ie=pG->ie;
  int j,js=pG->js,je=pG->je;
  int k,ks=pG->ks;
  
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->Phi[ks-k][j][i] = 3.0*pG->Phi[ks-k+1][j][i] - 3.0*pG->Phi[ks-k+2][j][i] + pG->Phi[ks-k+3][j][i];
#else  /* LINEAR EXTRAPOLATION */
        pG->Phi[ks-k][j][i] = 2.0*pG->Phi[ks-k+1][j][i] - pG->Phi[ks-k+2][j][i];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* OPEN (VACUUM) boundary conditions, Outer x3 boundary
 */
static void obc_fft_Phi_ox3(GridS *pG)
{
  int i,is=pG->is,ie=pG->ie;
  int j,js=pG->js,je=pG->je;
  int k,ke=pG->ke;
  
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifdef QUADRATIC_EXTRAPOLATION
        pG->Phi[ke+k][j][i] = 3.0*pG->Phi[ke+k-1][j][i] - 3.0*pG->Phi[ke+k-2][j][i] + pG->Phi[ke+k-3][j][i];
#else  /* LINEAR EXTRAPOLATION */
        pG->Phi[ke+k][j][i] = 2.0*pG->Phi[ke+k-1][j][i] - pG->Phi[ke+k-2][j][i];
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn static void ProlongateLater(GridS *pGrid)
 *  \brief PROLONGATION boundary conditions.  
 *
 * Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(GridS *pGrid)
{
  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~550 lines */
/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ix1(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_Phi_ix1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

/* Pack only Phi into send buffer */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=is+(nghost-1); i++){
        *(pSnd++) = pG->Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ox1(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_Phi_ox1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

/* Pack only Phi into send buffer */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie-(nghost-1); i<=ie; i++){
        *(pSnd++) = pG->Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ix2(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_Phi_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

/* Pack only Phi into send buffer */
  for (k=ks; k<=ke; k++){
    for (j=js; j<=js+(nghost-1); j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ox2(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_Phi_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

/* Pack only Phi into send buffer */

  for (k=ks; k<=ke; k++){
    for (j=je-(nghost-1); j<=je; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ix3(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_Phi_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

/* Pack only Phi into send buffer */

  for (k=ks; k<=ks+(nghost-1); k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_Phi_ox3(GridS *pG)
 *  \brief PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_Phi_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

/* Pack only Phi into send buffer */

  for (k=ke-(nghost-1); k<=ke; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        *(pSnd++) = pG->Phi[k][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ix1(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_Phi_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
//    for (j=js; j<=js; j++){
//MODIFIED BY HAO GONG
    for (j=js; j<=je; j++){
      for (i=is-nghost; i<=is-1; i++){
        pG->Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ox1(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_Phi_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie+1; i<=ie+nghost; i++){
        pG->Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ix2(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_Phi_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
    for (j=js-nghost; j<=js-1; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ox2(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_Phi_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks; k<=ke; k++){
    for (j=je+1; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ix3(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_Phi_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ks-nghost; k<=ks-1; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_Phi_ox3(GridS *pG)
 *  \brief UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_Phi_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ke = pG->ke;
  int i,j,k;
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

/* Manually unpack the data from the receive buffer */

  for (k=ke+1; k<=ke+nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->Phi[k][j][i] = *(pRcv++);
      }
    }
  }

  return;
}

#endif /* MPI_PARALLEL */

#endif /* SELF_GRAVITY */
