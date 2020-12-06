#ifndef _UTIL_H
#define _UTIL_H

struct point {
  double t;
  double y;
};

double euclidean_norm(int n, double *x);

double largest_vector_component(int n, double *x);

int pivot_row_index(int n, double **A, int from_index);

void interchange_pivot_row(int k, int pivot_index, double **A);

void initialize_vector(int n, double *b);

void initialize_matrix(int n, double **A);

double *allocate_vector(int n);

double **allocate_matrix(int n);

struct point **allocate_data_points(int n);

void free_matrix(int n, double **A);

void free_data_points(int n, struct point **data_points);

int read_size();

double *read_vector(int n);

double **read_matrix(int n);

struct point **read_lsp_input(int *n, int *m);

void print_vector(int n, double *x);

void print_matrix(int n, double **A);

void print_data_points(int n, struct point **data_points);

#endif
