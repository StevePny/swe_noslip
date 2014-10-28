#include <stdlib.h>
#include <math.h>
//#include <stdio.h>

#include "shallow_water.h"

void initialize_geometry(int xdim, int ydim, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order){

  int i, j, k, neighbor;
  unsigned int **array_2d;

  /*
   *  allocate temporary array holding coordinates
   */
  array_2d = calloc(xdim, sizeof(unsigned int*));
  for(i=0; i<xdim; i++){
    array_2d[i]=calloc(ydim, sizeof(unsigned int));
  }

  /*
   *  fill with coordinates, first fill interior 
   *
   *  5x5 example:   24 17 16 15 23 
   *                 18  6  7  8 14 
   *                 19  3  4  5 13 
   *                 20  0  1  2 12 
   *                 21  9 10 11 22
   */
  k=0;
  /*
   * fill interior
   */
  for(j=1; j<ydim-1; j++){
    for(i=1; i<xdim-1; i++){
      array_2d[i][j] = k;
      k++;
    }
  }

  /*
   * fill southern boundary (except corners)
   */
  for(i=1; i<xdim-1; i++){
    array_2d[i][0] = k;
    k++;
  }

  /*
   * fill eastern boundary (except corners)
   */
  for(j=1; j<ydim-1; j++){
    array_2d[xdim-1][j] = k;
    k++;
  }

  /*
   * fill northern boundary (except corners)
   */
  for(i=xdim-2; i>0; i--){
    array_2d[i][ydim-1] = k;
    k++;
  }

  /*
   * fill western boundary (except corners)
   */
  for(j=ydim-2; j>0; j--){
    array_2d[0][j] = k;
    k++;
  }

  /*
   * fill corners
   */
  array_2d[0][0] = k;
  k++;
  array_2d[xdim-1][0] = k;
  k++;
  array_2d[xdim-1][ydim-1] = k;
  k++;
  array_2d[0][ydim-1] = k;

  /*k=0;
  for(j=0; j<ydim; j++){
    for(i=0; i<xdim; i++){
      array_2d[i][j] = k;
      k++;
    }
  }*/
  /*
   *  fill the lat/lon arrays
   */
  k=0;
  for(j=0; j<ydim; j++){
    for(i=0; i<xdim; i++){
      lat[array_2d[i][j]] = j;
      lon[array_2d[i][j]] = i;
      print_out_order[k]=array_2d[i][j];
//    printf(" print order %i: %i\n", k, array_2d[i][j]);
      k++;
    }
  }

  /*
   * fill the neighbors array
   */
  for(k=0; k<xdim*ydim; k++){

    /*
     * get 2D coordinates
     */
    i=lon[k];
    j=lat[k];

    /*
     * determine neighbor in positive x-direction
     */
    neighbor=(i+1)%xdim;
    n[k][0] = array_2d[neighbor][j];

    /*
     * determine neighbor in positive y-direction
     */
    neighbor=(j+1)%ydim;
    n[k][1] = array_2d[i][neighbor];

    /*
     * determine neighbor in negative x-direction
     */
    neighbor=(i-1+xdim)%xdim;
    n[k][2] = array_2d[neighbor][j];

    /*
     * determine neighbor in negative y-direction
     */
    neighbor=(j-1+ydim)%ydim;
    n[k][3] = array_2d[i][neighbor];
  }

  //for(j=ydim-1; j>=0; j--){
  //  for(i=0; i<xdim; i++){
  //    printf("%u ", array_2d[i][j]);
  //  }
  //  printf("\n");
  //}

  //for(j=0; j<4; j++){
  //  for(i=0; i<xdim*ydim; i++){
  //    printf("%i ", n[i][j]);
  //  }
  //  printf("\n");
  //}

  //for(i=0; i<xdim*ydim; i++){
  //  printf("%u ", lat[i]);
  //}
  //printf("\n");

  //for(i=0; i<xdim*ydim; i++){
  //  printf("%u ", lon[i]);
  //}
  //printf("\n");

  //for(i=0; i<xdim*ydim; i++){
  //  printf("%u ", print_out_order[i]);
  //}
  //printf("\n");
  //exit(7);
}

/*----------------------------------------------------------------------------*/

void initialize_fields (double *u, double *v, double *P, double *x_forcing, double *y_forcing, double *z_forcing, int xdim, int ydim, double dx, double dy, double EquilibriumDepth, double A, int **n, unsigned int *lat, unsigned int *lon, enum bc_t bc)
{

  int i;
  int tdim;
  double eps;

  const double EL = ydim*dy;
  const double PI = 3.1415926;
  const double twoPi = PI+PI;
  const double DI = twoPi/(double)xdim;
  const double DJ = twoPi/(double)ydim;
  const double PCF = PI*PI*A*A/(EL*EL);
  double *PSI;     

  tdim = xdim*ydim;

  PSI = calloc(tdim, sizeof(double));

  for (i=0; i<tdim; i++) {

    PSI[i] = A*sin((lon[i]+0.5)*DI)*sin((lat[i]+0.5)*DJ);

    P[i] = PCF*(cos(4.0*((double)lon[i])*DI)+cos(4.0*((double)lat[i])*DJ))+EquilibriumDepth;
//    P[i] = PCF*(cos(2.0*lon[i]*DI)+cos(2.0*lat[i]*DJ))+EquilibriumDepth;
//    P[i] = 100.0* PCF*( sin(0.5*DI*lon[i]) + sin(0.5*DJ*lat[i])) + EquilibriumDepth;

//    P[i] =  lat[i] +EquilibriumDepth;
//    x_forcing[i] = cos(DJ*lat[i]);
    x_forcing[i] = sin(0.5*DJ*((double)lat[i]+0.5));
    x_forcing[i] = 0.0;
// I BELIEVE YOU MUST ALTER f_calc.c IF YOU INCLUDE Y or Z FORCING L(INE ~60 & 250)
    y_forcing[i] = 0.0;
    z_forcing[i] = 0.0;
  }

  /*
   *  initialize velocities
   */
  for (i=0;i<tdim;i++){
      u[i] = -(PSI[n[i][1]]-PSI[i])/dy;
      v[i] =  (PSI[n[i][0]]-PSI[i])/dx;
//        u[i] = 0.0;
//        v[i] = 0.0;
      //u[n[i][0]] = -(PSI[n[n[i][0]][1]]-PSI[n[i][0]])/dy;
      //v[n[i][1]] =  (PSI[n[n[i][0]][1]]-PSI[n[i][1]])/dx;
  }

  /*
   * handle no_slip boundary constditions
   */
  eps=1e-10;
  if(bc == no_slip){

    /*
     * add a small number in order to make the solver stable
     */
    //for(i=0; i<(xdim-2)*(ydim-2); i++){
    for(i=0; i<xdim*ydim; i++){
      u[i] += eps;
      v[i] += eps;
      P[i] += eps;
     }

    /*
     * set boundary velocities to zero for no_slip
     */
    for(i=(xdim-2)*(ydim-2); i<tdim; i++){
      u[i] =0;
      v[i] =0;
    }
  }
    free (PSI);
}
