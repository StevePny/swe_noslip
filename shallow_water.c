#include <stdlib.h>
#include <stdio.h>

#include "shallow_water.h"

int main (int argc, char ** argv)
{
  int i;
    
  /*
   *  this is set to either euler or runge_kutta_4
   */
  enum solver_t solver;

  /*
   *  bc (boundary condition) can be either periodic or no_slip
   */
  enum bc_t bc;

  /*
   * these specify the geometry of the area to be simulated 
   * tdim ist the total number of points, i.e. tdim = xdim *ydim 
   */
  int xdim, ydim, tdim;
  double dx, dy;

  /*
   * this holds the neighbors 
   */
  int **neighbors;

  /*
   * this holds latitude and longitude values of each point
   */
  unsigned int *lat, *lon;

  /*
   * this holds order of the data to be printed out
   */
  unsigned int *print_out_order;

  /*
   * these control temporal aspects of the simulation 
   */
  long int ncycle;
  double dt = 5.1;
  long int n_iter = 25800;
  int write_out_interval = 5000;

  /*
   * these hold the physical parameters 
   *
   * f0:              coriolis parameter (1/1sec)
   * beta:            linear beta for coriolis force (1/(meter sec))
   * forcingTerm:     Wind stress amplitude "tau_0"(N/meter^2)
   * dissipationTerm: "A" viscosity coefficient (meters^2/sec)
   * RayleighFriction: Rayleigh friction parameter (1/sec)
   */
  double f0 = 5.0e-5;
  double beta =2e-11;
  double forcingTerm = 0.0005;
  double dissipationTerm = 0.00005;
  double RayleighFriction = 5e-8;

  /*
   * upper layer equilibrium depth (meters)
   */
  double EquilibriumDepth = 50000.0;
  double A;

  /*
   * this holds the physical parameters and the forcing
   */
  double *parameters;
  double *x_forcing, *y_forcing, *z_forcing;

  /*
   * "fields" holds the values of the u, v and P fields
   */
  double *fields;
  double *u, *v, *P;     

  /*
   * "fields_prev" holds the u, v and P values of the previous time step
   */
  //double *fields_prev;
  //double *u_prev, *v_prev, *P_prev;     

  /*
   * these are temporary storage locations 
   * for the Runge-Kutta scheme
   */
  double *temp_dots_rk;
  double *temp_fields_rk;

  /*
   * read parameters from command line
   */
  if(argc == 2){

    get_parameters(argv[1], &xdim, &ydim, &dx, &dy, &n_iter, &dt, &write_out_interval, &bc, &solver, &f0, &beta, &forcingTerm, &dissipationTerm, &RayleighFriction, &EquilibriumDepth, &A);

  }
  else{

    fprintf(stderr, "ERROR: You need to give a file containing the parameters as an argument!\n");
    exit(1);

  }

  print_parameters(xdim, ydim, dx, dy, n_iter, dt, write_out_interval, bc, solver, f0, beta, forcingTerm, dissipationTerm, RayleighFriction, EquilibriumDepth, A);

  tdim = xdim*ydim;

  /*
   *  allocate arrays containing geometrical information and fill these
   */
  lat = calloc(tdim, sizeof(unsigned int)); 
  lon = calloc(tdim, sizeof(unsigned int)); 
  print_out_order = calloc(tdim, sizeof(unsigned int)); 
  neighbors = calloc(tdim, sizeof(int*)); 
  for(i=0; i< tdim; i++){
    neighbors[i] = calloc(4, sizeof(int));
  }

  initialize_geometry(xdim, ydim, neighbors, lat, lon, print_out_order);

  /*
   * allocate memory for physical parameters and forcing
   */
  parameters = calloc(5+3*tdim, sizeof(double));

  parameters[0] = f0;
  parameters[1] = beta;
  parameters[2] = forcingTerm;
  parameters[3] = dissipationTerm;
  parameters[4] = RayleighFriction;

  x_forcing = parameters + 5;
  y_forcing = parameters + 5 + tdim;
  z_forcing = parameters + 5 + 2*tdim;

  /*
   * allocate "fields" and set addresses for u, v, and P
   */
  fields = calloc(3*tdim, sizeof(double));
  u = fields;
  v = fields + tdim;
  P = fields + 2 * tdim;

  /*
   * allocate "fields_prev" and set addresses for u_prev, v_prev, and P_prev
   */
  //fields_prev = calloc(3*tdim, sizeof(double));
  //u_prev = fields_prev;
  //v_prev = fields_prev + tdim;
  //P_prev = fields_prev + 2 * tdim;

  double *fields_dot;
  fields_dot = calloc(3*tdim, sizeof(double));

  /*
   * allocate memory for temporary Runge-Kutta storage locations 
   */
  temp_dots_rk = calloc(3*tdim, sizeof(double));
  temp_fields_rk = calloc(3*tdim, sizeof(double));

  /*
   * initialize fields
   */
  initialize_fields (u, v, P, x_forcing, y_forcing, z_forcing, xdim, ydim, dx, dy, EquilibriumDepth, A, neighbors, lat, lon, bc);

  /*
   *  loop through time steps to do the actual simulation
   */
  for (ncycle=0; ncycle<n_iter; ncycle++) {
    //printf("step #%li\n", ncycle);

    /*
     * print result
     */
    if(ncycle % write_out_interval == 0){
      print_field(u, "u", ncycle, xdim, ydim, print_out_order);
      print_field(v, "v", ncycle, xdim, ydim, print_out_order);
      print_field(P, "P", ncycle, xdim, ydim, print_out_order);
      printf("%li ", ncycle);
      fflush(stdout);
    }
    Fcalc(fields_dot, fields, parameters, xdim, ydim, dx, dy, neighbors, lat, bc, print_out_order, ncycle);
 

    //for (i=0;i<3*tdim;i++) {
    //  fields_prev[i] = fields[i];
    //}

    if(solver == runge_kutta_4){

      /*
       *  Runge-Kutta 4th order
       *
       *                 dt
       *    y_1 = y_0 + ---- (y'_0 + 2 * y'_A + 2 * y'_B +y'_C)
       *                  6
       */

      /*                          dt
       * first step: y_A = y_0 + ---- y'_0
       *                          2
       */
      for (i=0; i < 3*tdim;  i++) {
        temp_fields_rk[i] = fields[i]+0.5*dt*fields_dot[i];
      }
      Fcalc(temp_dots_rk, temp_fields_rk, parameters, xdim, ydim, dx, dy, neighbors, lat, bc, print_out_order, ncycle);
      for (i=0; i < 3*tdim;  i++) {
        fields_dot[i] += 2*temp_dots_rk[i];
      }

      /*                           dt
       * second step: y_B = y_0 + ---- y'_A
       *                            2
       */
      for (i=0; i < 3*tdim;  i++) {
        temp_fields_rk[i] = fields[i]+0.5*dt*temp_dots_rk[i];
      }
      Fcalc(temp_dots_rk, temp_fields_rk, parameters, xdim, ydim, dx, dy, neighbors, lat, bc, print_out_order, ncycle);
      for (i=0; i < 3*tdim;  i++) {
        fields_dot[i] += 2*temp_dots_rk[i];
      }

      /*
       * third step: y_C = y_0 + dt * y'_B
       */
      for (i=0; i < 3*tdim;  i++){
        temp_fields_rk[i] = fields[i]+dt*temp_dots_rk[i];
      }
      Fcalc(temp_dots_rk, temp_fields_rk, parameters, xdim, ydim, dx, dy, neighbors, lat, bc, print_out_order, ncycle);
      for (i=0; i < 3*tdim;  i++) {
        fields_dot[i] += temp_dots_rk[i];
      }

      /*                          dt
       * final step: y_1 = y_0 + ---- (y'_0 + 2 * y'_A + 2 * y'_B +y'_C)
       *                           6
       */
      for (i=0; i < 3*tdim;  i++) {
        fields_dot[i] /= 6.0;
        fields[i] += dt*fields_dot[i];
      }
    }
    else{

      if(solver == euler){

        /*
         *  the Euler scheme
         */
        for (i=0; i < 3*tdim;  i++) {
          fields[i] += dt*fields_dot[i];
        }
      }
      else{
        fprintf(stderr, "ERROR: No valid solver was found\n");
        exit(2);
      }
    }
  }
  printf("\n");

  return(0);
}
