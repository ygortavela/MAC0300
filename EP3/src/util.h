#ifndef _UTIL_H
#define _UTIL_H

#define EPSILON 1e-6

struct point {
  double t;
  double y;
};

int backrow(int n, double **A, double *b);

double euclidean_norm(int n, int init, double *x);

double euclidean_norm_with_scaling(int n, int init, double *x);

double largest_vector_component(int n, int init, double *x);

int pivot_row_index(int n, double **A, int from_index);

void interchange_pivot_row(int k, int pivot_index, double **A);

void initialize_vector(int n, double *b);

void initialize_matrix(int n, int m, double **A);

double *allocate_vector(int n);

double **allocate_matrix(int n, int m);

struct point **allocate_data_points(int n);

void free_matrix(int n, double **matrix);

void free_data_points(int n, struct point **data_points);

int read_size();

struct point **read_lsp_input(int *n, int *m);

void print_vector(int n, double *x);

void print_matrix(int n, int m, double **matrix);

void print_data_points(int n, struct point **data_points);

double **build_transpose_coefficient_matrix_on_std_basis(int n, int m, struct point **data_points);

#endif
