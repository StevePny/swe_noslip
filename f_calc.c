#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "shallow_water.h"

/*
 * the variable names follow
 * R. Sadourny, "The Dynamics of Finite-Difference Models of the
 * Shallow-Water Equations", Journal of the Atmospheric Sciences,
 * Vol 32, p. 680, (1975)
 */

void Fcalc(double *fields_dot, const double *fields, const double *parameters, int xdim, int ydim, double dx, double dy, int **n, unsigned int *lat, enum bc_t bc, unsigned int *print_out_order, long int ncycle)
{

  int i, i_max, tdim;
  int southern_edge, eastern_edge, northern_edge, western_edge, sw_corner;

  tdim = xdim*ydim;
  if(bc == periodic){
    i_max = xdim*ydim;
  }
  if(bc == no_slip){
    i_max = (xdim-2)*(ydim-2);
  }
  southern_edge = (xdim-2)*(ydim-2);
  eastern_edge  = southern_edge+xdim-2;
  northern_edge = eastern_edge+ydim-2;
  western_edge  = northern_edge+xdim-2;
  sw_corner     = western_edge+ydim-2;
  //printf("FGD  %u %i %i\n", bc, tdim, i_max);
  //printf("FGE  S%i E%i N%i W%i\n", southern_edge, eastern_edge, northern_edge, western_edge);

  /*
  * some recurring numbers
   */
  const double S8 = 1.0/8.0;
  const double SX = 1.0/dx;
  const double SY = 1.0/dy;
  const double FSDX = 4.0/dx;
  const double FSDY = 4.0/dy;

  /*
   * these hold the physical parameters 
   *
   * f0:              coriolis parameter (1/1sec)
   * beta:            linear beta for coriolis force (1/(meter sec))
   * forcingTerm:     Wind stress amplitude "tau_0"(N/meter^2)
   * dissipationTerm: "A" viscosity coefficient (meters^2/sec)
   * RayleighFriction: Rayleigh friction parameter (1/sec)
   */
  double f0 = parameters[0];
  double beta = parameters[1];
  double forcingTerm = parameters[2];
  double dissipationTerm = parameters[3];
  double RayleighFriction = parameters[4];

  const double *forcingShapeX = parameters + 5;
  //const double *forcingShapeY = parameters + 5 + tdim;
  //const double *forcingShapeZ = parameters + 5 + 2*tdim;

  double *U, *V, *H, *eta;

  const double *u = fields;
  const double *v = fields + tdim;
  const double *P = fields + 2*tdim;

  double *u_dot = fields_dot;
  double *v_dot = fields_dot +tdim;
  double *P_dot = fields_dot + 2*tdim;

  /*
   * allocate memory for U, V, H and eta fields and set them to zero
   */
  U   = calloc(tdim, sizeof(double));
  V   = calloc(tdim, sizeof(double));
  H   = calloc(tdim, sizeof(double));
  eta = calloc(tdim, sizeof(double));

  /*     
   * compute U, V, eta and H
   */     
  for (i=0;i<i_max;i++) {

    U[i] = .5*(P[i]+P[n[i][2]])*u[i];
    V[i] = .5*(P[i]+P[n[i][3]])*v[i];

    eta[i] = (FSDX*(v[i]-v[n[i][2]])-FSDY*(u[i]-u[n[i][3]]))
            /(P[n[n[i][2]][3]]+P[n[i][3]]+P[i]+P[n[i][2]]);

    H[i] = P[i]+.25*(u[n[i][0]]*u[n[i][0]]+u[i]*u[i]
                    +v[n[i][1]]*v[n[i][1]]+v[i]*v[i]);
  }
  //print_field(U, "U", ncycle, xdim, ydim, print_out_order);
  //print_field(V, "V", ncycle, xdim, ydim, print_out_order);
  //print_field(eta, "eta", ncycle, xdim, ydim, print_out_order);
  //print_field(H, "H", ncycle, xdim, ydim, print_out_order);

  /*
   * in case of no_slip boundary conditions determine H and eta 
   * on the southern and eastern boundary, U and V remain zero
   */
  if(bc == no_slip){

    /*
     * update H and eta on southern edge, set U and V to zero
     */
    for(i=southern_edge; i<eastern_edge; i++){
      eta[i] = (-FSDY*u[n[i][1]])
        /(P[n[n[i][2]][1]]+P[n[i][1]]+P[i]+P[n[i][2]]);
      H[i]   = P[i]+0.25*v[n[i][1]]*v[n[i][1]];
      //U[i] = 0;
      //V[i] = 0;
    }

    /*
     * update H and eta on eastern edge, set U and V to zero
     */
    for(i=eastern_edge; i<northern_edge; i++){
      eta[i] = -FSDX*v[n[i][2]]/(P[n[n[i][2]][3]]+P[n[i][3]]+P[i]+P[n[i][2]]);
      H[i]   = P[i]+0.25*u[n[i][2]]*u[n[i][2]];
      //U[i] = 0;
      //V[i] = 0;
    }

    /*
     * update H and eta on northern edge, set U and V to zero
     */
    for(i=northern_edge; i<western_edge; i++){
      eta[i] = FSDY*u[n[i][3]]/(P[n[n[i][2]][3]]+P[n[i][3]]+P[i]+P[n[i][2]]);
      H[i]   = P[i]+0.25*v[n[i][3]]*v[n[i][3]];
    }

    /*
     * update H and eta on western edge, set U and V to zero
     */
    for(i=western_edge; i<sw_corner; i++){
      eta[i] = FSDX*v[n[i][0]]/(P[n[n[i][0]][3]]+P[n[i][0]]+P[n[i][3]]+P[i]);
      H[i]   = P[i]+0.25*u[n[i][0]]*u[n[i][0]];
    }

    /*
     * set eta, H, U and V to zero on the corners, this is actually redundant
     */
    for(i=sw_corner; i<tdim; i++){
      eta[i] = 0.0;
      H[i]   = 0.0;
      U[i]   = 0.0;
      V[i]   = 0.0;
    }
  }
  //print_field(eta, "eta_b", ncycle, xdim, ydim, print_out_order);
  //print_field(H, "H_b", ncycle, xdim, ydim, print_out_order);

  /*
   * use computed U, V, eta and H to compute 
   * the vector fields u_dot, v_dot and P_dot
   */
  for (i=0; i<i_max; i++){

    //u_dot[i] = 
    //  S8*(eta[i]+eta[n[i][1]])*(V[i]+V[n[i][0]]+V[n[i][1]]+V[n[n[i][0]][1]])
    //  -SX*(H[n[i][0]]-H[i]);
    //u_dot[n[i][0]] = 
    //    S8*(eta[n[n[i][0]][1]]+eta[n[i][0]])*(V[n[n[i][0]][1]]+V[n[i][1]]
    //    +V[i]+V[n[i][0]])-SX*(H[n[i][0]]-H[i]);
    u_dot[i] = 
        S8*(eta[n[i][1]]+eta[i])*(V[n[i][1]]+V[n[n[i][2]][1]]
        +V[n[i][2]]+V[i])-SX*(H[i]-H[n[i][2]]);

    //v_dot[i] = 
    //  -S8*(eta[i]+eta[n[i][0]])*(U[i]+U[n[i][0]]+U[n[i][1]]+U[n[n[i][0]][1]])
    //  -SY*(H[n[i][1]]-H[i]);
    //v_dot[n[i][1]] = -S8*(eta[n[n[i][0]][1]]+eta[n[i][1]])
    //  *(U[n[n[i][0]][1]]+U[n[i][1]]+U[i]+U[n[i][0]])
    //  -SY*(H[n[i][1]]-H[i]);
    v_dot[i] = -S8*(eta[n[i][0]]+eta[i])
      *(U[n[i][0]]+U[i]+U[n[i][3]]+U[n[n[i][3]][0]])
      -SY*(H[i]-H[n[i][3]]);

    P_dot[i] = -SX*(U[n[i][0]]-U[i]) -SY*(V[n[i][1]]-V[i]);
    int iuz=0;
    if(isnan(u_dot[i])){
      printf("u blew up %i\n", i);
      printf("gtrf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", eta[i], eta[n[i][1]], V[i], V[n[i][0]], V[n[i][1]], V[n[n[i][0]][1]], H[n[i][0]], H[i]);
      iuz=1;
    }
    if(isnan(v_dot[i])){
      printf("v blew up %i\n", i);
      iuz=1;
    }
    if(isnan(P_dot[i])){
      printf("P blew up %i\n", i);
      printf("vads   %lf  %lf  %lf %lf\n", U[n[i][0]], U[i], V[n[i][1]], V[i]);
      iuz=1;
    }
    if(iuz == 1){
      exit(8);
    }
  }

  /*
   * include forcing and dissipation - include "del squared" dissipation 
   * on the quantities u, v and P
   */
  for (i=0;i<i_max; i++){

      /*
       * include the coriolis force 
       */
      u_dot[i] += (f0 + beta*dy*lat[i])*v[i];
      v_dot[i] -= (f0 + beta*dy*lat[i])*u[i];

      /*
       * include the lateral (viscous) dissipation
       */
      u_dot[i] += dissipationTerm*
        ( (SX*SX) * (u[n[i][0]]+u[n[i][2]]-2.0*u[i])
         +(SY*SY) * (u[n[i][1]]+u[n[i][3]]-2.0*u[i]));

      v_dot[i] += dissipationTerm*
        ( (SX*SX) * (v[n[i][0]]+v[n[i][2]]-2.0*v[i])
         +(SY*SY) * (v[n[i][1]]+v[n[i][3]]-2.0*v[i]));

      /*
       * bottom friction
       */
      u_dot[i] -= RayleighFriction*u[i];
      v_dot[i] -= RayleighFriction*v[i];

      /*
       * and forcing
       */
      u_dot[i] += forcingTerm*forcingShapeX[i];
      //v_dot[i] += forcingTerm*forcingShapeY[i];
      //P_dot[i] += forcingTerm*forcingShapeZ[i];
  }
  free(U);
  free(V);
  free(H);
  free(eta);
}
