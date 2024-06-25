/*--- Reference: https://nbviewer.org/github/barbagroup/CFDPython/blob/master/lessons/05_Step_4.ipynb ---*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct DOMAIN_s {
  int nx;
  int nt;
  double dx;
  double dt;
  double *x;
} domain_t;

typedef struct FIELD_s {
  double nu;
  double sigma;
  double *u;
  double *u_prev;
} field_t;

int copy_field(double *src, double *dest, int n)
{
  int i;
  for(i = 0; i < n; ++i) {
    dest[i] = src[i];
  }

  return 0;
}

void update_prev_field(field_t *field, domain_t *dom) 
{
  copy_field(field->u, field->u_prev, dom->nx);
}


domain_t* create_domain() 
{
  domain_t *dom;

  dom = (domain_t*) malloc(sizeof(domain_t));

  dom->nx = 0;
  dom->dx = 0.0;
  dom->nt = 0;
  dom->dt = 0.0;
  dom->x = NULL;

  return dom;
}

int setup_domain_spatial_points(domain_t *dom) 
{
  int i;

  if(dom->x != NULL) {
    fprintf(stdout, "[ERROR] Spatial points array not empty!\n");
    return -1;
  }

  if(!(dom->x = (double*) malloc(sizeof(double) * dom->nx))) {
    return -1;
  }

  for(i = 0; i < dom->nx; ++i) {
    dom->x[i] = i * dom->dx;
  }
   
  return 0;
}

int set_domain(domain_t *dom, int nx, double dx, int nt, double dt)
{
  dom->nx = nx;
  dom->dx = dx;
  dom->nt = nt;
  dom->dt = dt;

  setup_domain_spatial_points(dom);

  return 0;
}

double get_phi(double x, double t, double nu)
{
  return exp(-(x - 4*t) * (x - 4*t) / (4*nu * (t + 1))) + \
         exp(-(x - 4*t -2*M_PI) * (x - 4*t -2*M_PI) / (4*nu*(t + 1)));
}

double get_phi_deriv(double x, double t, double nu)
{
  return -(2*x - 8*t) * \
         exp(-(x - 4*t)*(x - 4*t) / (4*nu*(t + 1))) / \
         (4*nu*(t + 1)) - \
         (2*x - 8*t - 4*M_PI) * \
         exp(-(x - 4*t - 2*M_PI) * (x - 4*t - 2*M_PI) / \
         (4*nu*(t + 1))) / (4*nu*(t + 1));
}

double get_vel(double x, double t, double nu)
{
  return ((-2 * nu * get_phi_deriv(x, t, nu)) / get_phi(x, t, nu)) + 4;
}

int set_field_initial_conditions(domain_t *dom, field_t *field)
{
  int i;
  for(i = 0; i < dom->nx; ++i) {
    field->u[i] = get_vel(dom->x[i], 0, field->nu);
  }
}

field_t* create_field(domain_t *dom, double nu, double sigma) 
{
  field_t *field;
  
  field = (field_t*) malloc(sizeof(field_t));
  if(!(field->u = (double*) calloc(dom->nx, sizeof(double)))) {
    return NULL;
  }

  if(!(field->u_prev = (double*) calloc(dom->nx, sizeof(double)))) {
    return NULL;
  }

  field->nu = nu;
  field->sigma = sigma;

  set_field_initial_conditions(dom, field);
  update_prev_field(field, dom);

  return field;
}

int free_domain(domain_t *dom) 
{
  free(dom->x);
  free(dom);

  return 0;
}

int free_field(field_t *field)
{ free(field->u);
  free(field);
}

int set_boundary_condition(field_t *field, domain_t *dom)
{
  int i, n;
  double c;
  double dt, dx, nu;
  double *u, *un;

  u = field->u;
  un = field->u_prev;
  n = dom->nx;
  dt = dom->dt;
  dx = dom->dx;
  nu = field->nu;
  
  u[n - 1] = un[n - 1] - 
         un[n - 1] * (dt/dx) * (un[n - 1] - un[n - 2]) + 
         nu * (dt/(dx * dx)) * (un[1] + un[n - 2] - 2 * un[n - 1]);
  u[0] = u[n - 1];
}

void execute_time_step(double (*f)(field_t*, int), domain_t *dom, field_t *field)
{
  int i, n;
  double c;
  double dt, dx, nu;
  double *u, *un;

  u = field->u;
  un = field->u_prev;
  n = dom->nx;
  dt = dom->dt;
  dx = dom->dx;
  nu = field->nu;

  for(i = 1; i < n - 1; ++i) {
    c = (double)(*f)(field, i);
    u[i] = un[i] - 
           un[i] * (dt/dx) * (un[i] - un[i - 1]) + 
           nu * (dt/(dx * dx)) * (un[i + 1] + un[i - 1] - 2 * un[i]);
    set_boundary_condition(field, dom);
  }
  update_prev_field(field, dom);
}

double advection_velocity(field_t *field, int i) 
{
  return 1.0;
  //return field->u_prev[i]; 
}

void write_to_file(char* file_name, double *arr, int n) 
{
  int i;
  FILE *fp;

  fp = fopen(file_name, "a");
  if (fp == NULL) {
    fprintf(stdout, "[ERROR] Failed to open file: %s\n", file_name);
    return;
  }
  for (i = 0; i < n; ++i) {
    fprintf(fp, "%e,", arr[i]);
  }
  fprintf(fp, "\n");
  fclose(fp);
}

void remove_file_if_exists(char* file_name)
{
  if(remove(file_name) == 0) {
    fprintf(stdout, "Output file removed successfully!\n");
  }
  else {
    fprintf(stdout, "[ERROR] Failed to delete output file...\n");
  }

}

int set_domain_dt_from_field(domain_t *dom, field_t *field)
{
  dom->dt = field->sigma * dom->dx * field->nu;
  return 0;
}

int main()
{
  domain_t *dom;
  field_t *field;
  int i, num_iters;
  char *file_name;

  file_name = "step4_out.csv";
  int n_sp_pts = 101;

  dom = create_domain();
  set_domain(dom, n_sp_pts, (2.0 * M_PI)/(n_sp_pts - 1), 250, -1);
  field = create_field(dom, 0.07, 1.0);
  set_domain_dt_from_field(dom, field);

  remove_file_if_exists(file_name);
  for(i = 0; i < dom->nt; ++i) {
    write_to_file(file_name, field->u, dom->nx);
    execute_time_step(advection_velocity, dom, field);
  }

  free_field(field);
  free_domain(dom);
  return 0;
}
