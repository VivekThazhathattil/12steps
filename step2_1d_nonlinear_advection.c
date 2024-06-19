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

field_t* create_field(domain_t *dom) 
{
  int i;
  field_t *field;
  
  field = (field_t*) malloc(sizeof(field_t));
  if(!(field->u = (double*) calloc(dom->nx, sizeof(double)))) {
    return NULL;
  }

  if(!(field->u_prev = (double*) calloc(dom->nx, sizeof(double)))) {
    return NULL;
  }

  for(i = 0; i < dom->nx; ++i) {
    if(dom->x[i] <= 1.0 && dom->x[i] >= 0.5) {
      field->u[i] = 3.0;
    }
    else {
      field->u[i] = 1.0;
    }
  }

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

void execute_time_step(double (*c)(field_t*, int), domain_t *dom, field_t *field)
{
  int i;
  double vel;
  for(i = 1; i < dom->nx; ++i) {
    vel = (double)(*c)(field, i);
    field->u[i] = field->u_prev[i] - \
       vel * (dom->dt/dom->dx) * (field->u_prev[i] - field->u_prev[i - 1]);
  }
  update_prev_field(field, dom);
}

double advection_velocity(field_t *field, int i) 
{
  return field->u_prev[i]; 
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

int main()
{
  domain_t *dom;
  field_t *field;
  int i, num_iters;
  char *file_name;

  num_iters = 500;
  file_name = "step2_out.csv";

  dom = create_domain();
  set_domain(dom, 101, 2.0/50, 20, 2.5e-3);
  field = create_field(dom);

  remove_file_if_exists(file_name);
  for(i = 0; i < num_iters; ++i) {
    write_to_file(file_name, field->u, dom->nx);
    execute_time_step(advection_velocity, dom, field);
  }

  free_field(field);
  free_domain(dom);
  return 0;
}
