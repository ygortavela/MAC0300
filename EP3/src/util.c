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

double euclidean_norm(int n, int init, double *x) {
  double dot_sum = 0;

  for (int i = init; i < n; i++) {
    dot_sum += x[i] * x[i];
  }

  return sqrt(dot_sum);
}

double euclidean_norm_with_scaling(int n, int init, double *x) {
  double normalized_component, normalized_dot_sum = 0, largest_component = largest_vector_component(n, init, x);

  for (int i = init; i < n; i++) {
    normalized_component = x[i]/largest_component;
    normalized_dot_sum += normalized_component * normalized_component;
  }

  return largest_component * sqrt(normalized_dot_sum);
}

double largest_vector_component(int n, int init, double *x) {
  double temp, max = fabs(x[init]);

  for (int i = init + 1; i < n; i++) {
    temp = fabs(x[i]);

    if (temp > max)
      max = temp;
  }

  return max;
}

int pivot_row_index(int n, int m, double **A, int init) {
  double largest_norm = euclidean_norm_with_scaling(n, init, A[init]), aux;
  int max_index = init;

  for (int i = init + 1; i < m; i++) {
    aux = euclidean_norm_with_scaling(n, init, A[i]);

    if (aux > largest_norm) {
      max_index = i;
      largest_norm = aux;
    }
  }

  return max_index;
}

void interchange_pivot_row(int k, int pivot_index, double **A) {
  double *temp = A[k];

  A[k] = A[pivot_index];
  A[pivot_index] = temp;
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
