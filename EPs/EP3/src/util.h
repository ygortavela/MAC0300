#ifndef _UTIL_H
#define _UTIL_H

#define EPSILON 1e-6

struct point {
  double t;
  double y;
};

int backrow(int n, double **A, double *b);

double squared_sum_of_vector(int n, int init, double *x);

double euclidean_norm(int n, int init, double *x);

int largest_vector_component_index(int n, int init, double *x);

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

double largest_matrix_component(int n, int m, double **A);

void system_rescale(int n, int m, double **A, struct point **data_points);

double *build_cached_row_norms_vector(int n, int m, double **A);

double frobenius_norm_using_cached_norms_vector(int n, double *cached_norms);

void interchange_pivot_row(int k, int pivot_index, double **A);

void update_cached_norms_vector(int n, int init, double *cached_norms, double **A);

void swap_int_vector_value(int k, int pivot_index, int *v);

void swap_double_vector_value(int k, int pivot_index, double *v);

void assemble_permuted_solution(int r, int n, double *x, int *permutation);

void print_polynomial_solution(int n, double *x);

#endif
