enum solver_t {euler, runge_kutta_4, undef_solver};
enum bc_t {periodic, no_slip, undef_bc};

/*
 * this is from f_calc.c
 */
void Fcalc(double *fields_dot, const double *fields, const double *parameters, int xdim, int ydim, double dx, double dy, int **n, unsigned int *lat, enum bc_t bc, unsigned int *print_out_order, long int ncycle);

/*
 * these are from init.c
 */
void initialize_geometry(int xdim, int ydim, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order);
void initialize_fields (double *u, double *v, double *P, double *x_forcing, double *y_forcing, double *z_forcing, int xdim, int ydim, double dx, double dy, double EquilibriumDepth, double A, int **n, unsigned int *lat, unsigned int *lon, enum bc_t bc);

/*
 * these are from io.c
 */
void get_parameters(char *file_name, int *xdim, int *ydim, double *dx, double *dy, long int *n_iter, double *dt, int *write_out_interval, enum bc_t *bc, enum solver_t *solver, double *f0, double *beta, double *forcingTerm, double *dissipationTerm, double *RayleighFriction, double *EquilibriumDepth, double *A);
void print_parameters(int xdim, int ydim, double dx, double dy, long int n_iter, double dt, int write_out_interval, enum bc_t bc, enum solver_t solver, double f0, double beta, double forcingTerm, double dissipationTerm, double RayleighFriction, double EquilibriumDepth, double A);
void print_field(double *field, char *name, int time_step, int xdim, int ydim, unsigned int *print_out_order);
