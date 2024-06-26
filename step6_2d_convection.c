#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct DOMAIN_s {
  int nx;
  int ny;
  double dx;
  double dy;
  double dt;
  int nt;
  double *x;
  double *y;
} domain_t;

typedef struct FIELD_s {
  double *u;
  double *u_prev;
  double *v;
  double *v_prev;
  double sigma;
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
  copy_field(field->u, field->u_prev, dom->nx * dom->ny);
  copy_field(field->v, field->v_prev, dom->nx * dom->ny);
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
  dom->y = NULL;

  return dom;
}

int setup_domain_spatial_points(domain_t *dom) 
{
  int i;

  if(dom->x != NULL || dom->y != NULL) {
    fprintf(stdout, "[ERROR] Spatial points array not empty!\n");
    return -1;
  }

  if(!(dom->x = (double*) malloc(sizeof(double) * dom->nx))) {
    return -1;
  }

  if(!(dom->y = (double*) malloc(sizeof(double) * dom->ny))) {
    return -1;
  }

  for(i = 0; i < dom->nx; ++i) {
    dom->x[i] = i * dom->dx;
  }

  for(i = 0; i < dom->ny; ++i) {
    dom->y[i] = i * dom->dy;
  }
   
  return 0;
}

int set_domain(domain_t *dom, int nx, int ny, double dx, double dy, int nt, double dt)
{
  dom->nx = nx;
  dom->ny = ny;
  dom->dx = dx;
  dom->dy = dy;
  dom->nt = nt;
  dom->dt = dt;

  setup_domain_spatial_points(dom);

  return 0;
}

int set_dt(domain_t *dom, field_t *field)
{
  dom->dt = field->sigma * dom->dx;
  return 0;
}

field_t* create_field(domain_t *dom) 
{
  int i, j;
  int nx, ny;
  field_t *field;
  
  nx = dom->nx;
  ny = dom->ny;

  field = (field_t*) malloc(sizeof(field_t));
  if(!(field->u = (double*) calloc(nx * ny, sizeof(double)))) {
    return NULL;
  }

  if(!(field->u_prev = (double*) calloc(nx * ny, sizeof(double)))) {
    return NULL;
  }

  if(!(field->v = (double*) calloc(nx * ny, sizeof(double)))) {
    return NULL;
  }

  if(!(field->v_prev = (double*) calloc(nx * ny, sizeof(double)))) {
    return NULL;
  }

  field->sigma = 0.2;

  for(i = 0; i < nx; ++i) {
    for(j = 0; j < ny; ++j) {
      if(dom->x[i] <= 1.0 && dom->x[i] >= 0.5 && \
         dom->y[j] <= 1.0 && dom->y[j] >= 0.5) 
      {
        field->u[i*ny + j] = 2.0;
        field->v[i*ny + j] = 2.0;
      }
      else 
      {
        field->u[i*ny + j] = 1.0;
        field->v[i*ny + j] = 1.0;
      }
    }
  }

  update_prev_field(field, dom);

  return field;
}

int free_domain(domain_t *dom) 
{
  free(dom->x);
  free(dom->y);
  free(dom);

  return 0;
}

int free_field(field_t *field)
{
  free(field->u);
  free(field);
}

void execute_time_step(double (*c)(int, field_t*, int), domain_t *dom, field_t *field)
{
  int i, j;
  double *un, *u;
  double *vn, *v;
  double cx, cy;
  double dt, dx, dy;
  int nx, ny;

  u = field->u;
  un = field->u_prev;
  v = field->v;
  vn = field->v_prev;

  dt = dom->dt;
  dx = dom->dx;
  dy = dom->dy;
  nx = dom->nx;
  ny = dom->ny;

  for(i = 1; i < nx; ++i) {
    for(j = 1; j < ny; ++j) {
      cx = (double)(*c)(1, field, i*ny + j);
      cy = (double)(*c)(2, field, i*ny + j);
      u[i*ny + j] = un[i*ny + j] \
              - cx * (dt/dx) * (un[i*ny + j] - un[(i - 1)*ny + j]) \
              - cy * (dt/dy) * (un[i*ny + j] - un[i*ny + (j - 1)]);
      v[i*ny + j] = vn[i*ny + j] \
              - cx * (dt/dx) * (vn[i*ny + j] - vn[(i - 1)*ny + j]) \
              - cy * (dt/dy) * (vn[i*ny + j] - vn[i*ny + (j - 1)]);
    }
  }

  update_prev_field(field, dom);
}

int set_boundary_conditions(domain_t *dom, field_t *field)
{
  int i, j;
  int nx, ny;
  double *u, *v;

  nx = dom->nx;
  ny = dom->ny;
  u = field->u;
  v = field->v;

  for(i = 0; i < nx; ++i) {
    u[i*ny + (0)] = 1.0;
    u[i*ny + (ny - 1)] = 1.0;
    v[i*ny + (0)] = 1.0;
    v[i*ny + (ny - 1)] = 1.0;
  }
  for(j = 0; j < ny; ++j) {
    u[(0*ny) + j] = 1.0;
    u[(nx - 1)*ny + j] = 1.0;
    v[(0*ny) + j] = 1.0;
    v[(nx - 1)*ny + j] = 1.0;
  }

  return 0;
}

double advection_velocity(int ndx, field_t* field, int i) 
{
  if(ndx == 1) {
    return field->u_prev[i];
  }
  return field->v_prev[i];
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
    fprintf(stdout, "File removed successfully!\n");
  }
  else {
    fprintf(stdout, "[ERROR] Failed to delete output file...\n");
  }

}

int main()
{
  domain_t *dom;
  field_t *field;
  int i, num_iters;
  char *file_name;

  file_name = "step6_out.csv";

  dom = create_domain();
  /*-- set_domain: dom, nx, ny, dx, dy, nt, dt --*/
  set_domain(dom, 81, 81, 2.0/80, 2.0/80, 200, -1);
  field = create_field(dom);
  set_dt(dom, field);

  remove_file_if_exists(file_name);
  for(i = 0; i < dom->nt; ++i) {
    write_to_file(file_name, field->u, dom->nx * dom->ny);
    execute_time_step(advection_velocity, dom, field);
    set_boundary_conditions(dom, field);
  }

  free_field(field);
  free_domain(dom);
  return 0;
}
