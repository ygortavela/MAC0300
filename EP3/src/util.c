#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "util.h"

int backrow(int n, double **A, double *b) {
  for (int j = n - 1; j >= 0; j--) {
    if (fabs(A[j][j]) < EPSILON) return -1;

    b[j] /= A[j][j];

    for (int i = j - 1; i >= 0; i--)
      b[i] -= A[j][i] * b[j];
  }

  return 0;
}

double squared_sum_of_vector(int n, int init, double *x) {
  double dot_sum = 0;

  for (int i = init; i < n; i++) {
    dot_sum += x[i] * x[i];
  }

  return dot_sum;
}

double euclidean_norm(int n, int init, double *x) {
  double dot_sum = squared_sum_of_vector(n, init, x);

  return sqrt(dot_sum);
}

double euclidean_norm_with_scaling(int n, int init, double *x) {
  int max_index = largest_vector_component_index(n, init, x);
  double normalized_component, normalized_dot_sum = 0, largest_component = x[max_index];

  for (int i = init; i < n; i++) {
    normalized_component = x[i]/largest_component;
    normalized_dot_sum += normalized_component * normalized_component;
  }

  return largest_component * sqrt(normalized_dot_sum);
}

int largest_vector_component_index(int n, int init, double *x) {
  double temp, max = fabs(x[init]);
  int index = init;

  for (int i = init + 1; i < n; i++) {
    temp = fabs(x[i]);

    if (temp > max) {
      max = temp;
      index = i;
    }
  }

  return index;
}

int pivot_row_index(int n_rows, int m_columns, double **A, int init) {
  double largest_norm = euclidean_norm_with_scaling(m_columns, init, A[init]), aux;
  int max_index = init;

  for (int i = init + 1; i < n_rows; i++) {
    aux = euclidean_norm_with_scaling(m_columns, init, A[i]);

    if (aux > largest_norm) {
      max_index = i;
      largest_norm = aux;
    }
  }

  return max_index;
}

void initialize_vector(int n, double *x) {
  for (int i = 0; i < n; i++)
    x[i] = 0;
}

void initialize_matrix(int n, int m, double **X) {
  for (int i = 0; i < n; i++)
    initialize_vector(m, X[i]);
}

double *allocate_vector(int n) {
  return malloc(n * sizeof(double));
}

double **allocate_matrix(int n, int m) {
  double **A = malloc(n * sizeof(double *));

  for (int i = 0; i < n; i++)
    A[i] = allocate_vector(m);

  return A;
}

struct point **allocate_data_points(int n) {
  struct point **data_points = malloc(n * sizeof(struct point *));

  for (int i = 0; i < n; i++)
    data_points[i] = malloc(sizeof(struct point));

  return data_points;
}

void free_matrix(int n, double **matrix) {
  for (int i = 0; i < n; i++)
    free(matrix[i]);

  free(matrix);
}

void free_data_points(int n, struct point **data_points) {
  for (int i = 0; i < n; i++)
    free(data_points[i]);

  free(data_points);
}

int read_size() {
  int n;

  scanf("%d", &n);

  return n;
}

struct point **read_lsp_input(int *n, int *m) {
  struct point **data_points;

  *m = read_size();
  // incrementing m to avoid changes on loops throught the code, remembering the fact that m is the degree of the polynomial
  *m += 1;
  *n = read_size();
  data_points = allocate_data_points(*n);

  for (int i = 0; i < *n; i++)
    scanf("%lf %lf", &data_points[i]->t, &data_points[i]->y);

  return data_points;
}

void print_vector(int n, double *x) {
  for (int i = 0; i < n; i++)
    printf("%lf ", x[i]);

  printf("\n");
}

void print_matrix(int n, int m, double **matrix) {
  for (int i = 0; i < n; i++)
    print_vector(m, matrix[i]);
}

void print_data_points(int n, struct point **data_points) {
  for (int i = 0; i < n; i++)
    printf("t_%d y_%d: %lf %lf\n", i, i, data_points[i]->t, data_points[i]->y);
}

double **build_transpose_coefficient_matrix_on_std_basis(int n, int m, struct point **data_points) {
  double **matrix;

  matrix = allocate_matrix(m, n);

  for (int j = 0; j < n; j++)
      matrix[0][j] = 1;

  for (int i = 1; i < m; i++) {
    for (int j = 0; j < n; j++)
      matrix[i][j] = matrix[i - 1][j] * data_points[j]->t;
  }

  return matrix;
}

double largest_matrix_component(int n, int m, double **A) {
  int max_index = largest_vector_component_index(m, 0, A[0]), aux_index;
  double max = A[0][max_index], aux;

  for (int i = 1; i < n; i++) {
    aux_index = largest_vector_component_index(m, 0, A[i]);
    aux = A[i][max_index];

    if (aux > max) {
      max_index = aux_index;
      max = aux;
    }
  }

  return max;
}

void system_rescale(int n, int m, double **A, struct point **data_points) {
  double largest_component = largest_matrix_component(n, m, A);

  for (int i = 0; i < n; i++) {
    data_points[i]->y /= largest_component;

    for (int j = 0; j < m; j++)
      A[i][j] /= largest_component;
  }
}

double *build_cached_row_norms_vector(int n, int m, double **A) {
  double *cached_norms = allocate_vector(n);

  for (int i = 0; i < n; i++)
    cached_norms[i] = squared_sum_of_vector(m, 0, A[i]);

  return cached_norms;
}

double frobenius_norm_using_cached_norms_vector(int n, double *cached_norms) {
  double sum = 0;

  for (int i = 0; i < n; i++)
    sum += cached_norms[i];

  return sqrt(sum);
}

void interchange_pivot_row(int k, int pivot_index, double **A) {
  double *temp = A[k];

  A[k] = A[pivot_index];
  A[pivot_index] = temp;
}

void update_cached_norms_vector(int n, int init, double *cached_norms, double **A) {
  if (init == 0) return;

  for (int i = init; i < n; i++)
    cached_norms[i] -= A[init - 1][i] * A[init - 1][i];
}

void interchange_cached_norms_values(int k, int pivot_index, double *cached_norms) {
  double temp = cached_norms[k];

  cached_norms[k] = cached_norms[pivot_index];
  cached_norms[pivot_index] = temp;
}
