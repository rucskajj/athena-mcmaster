#include "copyright.h"
/*==============================================================================
 * FILE: output_vtk.c
 *
 * PURPOSE: Function to write a single variable in VTK "legacy" format.  With
 *   SMR, dumps are made for all levels and domains, unless nlevel and ndomain
 *   are specified in <output> block.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_vtk() - writes VTK file (single variable).
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   output_vtk_2d() - write vtk file for 2D data
 *   output_vtk_3d() - write vtk file for 3D data
 *============================================================================*/

static void output_vtk_2d(MeshS *pM, OutputS *pOut, int nl, int nd);
static void output_vtk_3d(MeshS *pM, OutputS *pOut, int nl, int nd);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* output_vtk:   */

void output_vtk(MeshS *pM, OutputS *pOut)
{
  int nl,nd;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){

        if (pOut->ndim == 3) {
          output_vtk_3d(pM, pOut, nl, nd);
        } else if (pOut->ndim == 2) {
          output_vtk_2d(pM, pOut, nl, nd);
        } else {
          ath_error("[output_vtk]: Only able to output 2D or 3D");
        }
      }}
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* output_vtk_2d: writes 2D data  */

static void output_vtk_2d(MeshS *pM, OutputS *pOut, int nl, int nd)
{
  GridS *pGrid=pM->Domain[nl][nd].Grid;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int big_end = ath_big_endian();
  int nx1,nx2,ndata0;
  Real **data2d=NULL;
  double x1, x2, x3, dx1, dx2, dx3;

/* Allocate memory for and compute 2D array of data */
  data2d = OutData2(pGrid,pOut,&nx1,&nx2);
  if (data2d == NULL) return; /* data not in range of Grid */

/* construct output filename.  pOut->id will either be name of variable,
 * if 'id=...' was included in <ouput> block, or 'outN' where N is number of
 * <output> block.  */
  if (nl>0) {
    plev = &levstr[0];
    sprintf(plev,"lev%d",nl);
  }
  if (nd>0) {
    pdom = &domstr[0];
    sprintf(pdom,"dom%d",nd);
  }

  if((fname = ath_fname(plev,pM->outfilename,plev,pdom,num_digit,
      pOut->num,pOut->id,"vtk")) == NULL){
    ath_error("[output_vtk]: Error constructing filename\n");
  }

/* open output file */
  if((pfile = fopen(fname,"w")) == NULL){
    ath_error("[output_vtk]: Unable to open vtk file %s\n",fname);
  }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

  fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */

  fprintf(pfile,"Really cool Athena data at time= %e, level= %i, domain= %i\n",
    pGrid->time,nl,nd);

/*  3. File format */

  fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */
  x1 = pGrid->MinX[0];
  x2 = pGrid->MinX[1];
  x3 = pGrid->MinX[2];
  ndata0 = nx1*nx2;
  dx1 = (pOut->reduce_x1 == 1 ? pGrid->dx1 * pGrid->Nx[0] : pGrid->dx1);
  dx2 = (pOut->reduce_x2 == 1 ? pGrid->dx2 * pGrid->Nx[1] : pGrid->dx2);
  dx3 = (pOut->reduce_x3 == 1 ? pGrid->dx3 * pGrid->Nx[2] : pGrid->dx3);

  fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
  fprintf(pfile,"DIMENSIONS %d %d %d\n",nx1+1,nx2+1,1);
  fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
  fprintf(pfile,"SPACING %e %e %e \n",dx1,dx2,dx3);

/*  5. Data  */

  fprintf(pfile,"CELL_DATA %d \n", ndata0);

/* Write data */

  fprintf(pfile,"SCALARS %s float\n", pOut->id);
  fprintf(pfile,"LOOKUP_TABLE default\n");
  if(!big_end) ath_bswap(data2d[0],sizeof(float),ndata0);
  fwrite(data2d[0],sizeof(float),(size_t)ndata0,pfile);

/* close file and free memory */

  fclose(pfile);
  free_2d_array(data2d);
  return;
}

/*----------------------------------------------------------------------------*/
/* output_vtk_3d: writes 3D data  */

static void output_vtk_3d(MeshS *pM, OutputS *pOut, int nl, int nd)
{
  GridS *pGrid=pM->Domain[nl][nd].Grid;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int big_end = ath_big_endian();
  int nx1,nx2,nx3,ndata0, k;
  Real ***data3d=NULL;
  double x1, x2, x3;

/* construct output filename.  pOut->id will either be name of variable,
 * if 'id=...' was included in <ouput> block, or 'outN' where N is number of
 * <output> block.  */
  if (nl>0) {
    plev = &levstr[0];
    sprintf(plev,"lev%d",nl);
  }
  if (nd>0) {
    pdom = &domstr[0];
    sprintf(pdom,"dom%d",nd);
  }

  if((fname = ath_fname(plev,pM->outfilename,plev,pdom,num_digit,
      pOut->num,pOut->id,"vtk")) == NULL){
    ath_error("[output_vtk]: Error constructing filename\n");
  }

/* open output file */
  if((pfile = fopen(fname,"w")) == NULL){
    ath_error("[output_vtk]: Unable to open vtk file %s\n",fname);
  }

/* Allocate memory for and compute 3D array of data values */
  data3d = OutData3(pGrid,pOut,&nx1,&nx2,&nx3);

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

  fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */

  fprintf(pfile,"Really cool Athena data at time= %e, level= %i, domain= %i\n",
    pGrid->time,nl,nd);

/*  3. File format */

  fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */
  x1 = pGrid->MinX[0];
  x2 = pGrid->MinX[1];
  x3 = pGrid->MinX[2];
  ndata0 = nx1*nx2;

  fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
  fprintf(pfile,"DIMENSIONS %d %d %d\n",nx1+1,nx2+1,nx3+1);
  fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
  fprintf(pfile,"SPACING %e %e %e \n",pGrid->dx1,pGrid->dx2,pGrid->dx3);

/*  5. Data  */

  fprintf(pfile,"CELL_DATA %d \n", nx1*nx2*nx3);

/* Write data */

  fprintf(pfile,"SCALARS %s float\n", pOut->id);
  fprintf(pfile,"LOOKUP_TABLE default\n");
  for (k=0; k<nx3; k++) {
    if(!big_end) ath_bswap(data3d[k][0],sizeof(float),ndata0);
    fwrite(data3d[k][0],sizeof(float),(size_t)ndata0,pfile);
  }

/* close file and free memory */

  fclose(pfile);
  free_3d_array(data3d);
  return;
}
