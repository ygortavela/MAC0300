#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "util.h"

double euclidean_norm(int n, double *x) {
  double normalized_component, normalized_dot_sum = 0, largest_component = largest_vector_component(n, x);

  for (int i = 0; i < n; i++) {
    normalized_component = x[i]/largest_component;
    normalized_dot_sum += normalized_component * normalized_component;
  }

  return largest_component * sqrt(normalized_dot_sum);
}

double largest_vector_component(int n, double *x) {
  double temp, max = fabs(x[0]);

  for (int i = 1; i < n; i++) {
    temp = fabs(x[i]);

    if (temp > max)
      max = temp;
  }

  return max;
}

int pivot_row_index(int n, double **A, int from_index) {
  int max_index = from_index;

  for (int i = from_index + 1; i < n; i++) {
    if (A[i][from_index] > A[max_index][from_index])
      max_index = i;
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

void initialize_matrix(int n, double **X) {
  for (int i = 0; i < n; i++)
    initialize_vector(n, X[i]);
}

double *allocate_vector(int n) {
  return malloc(n * sizeof(double));
}

double **allocate_matrix(int n) {
  double **A = malloc(n * sizeof(double *));

  for (int i = 0; i < n; i++)
    A[i] = allocate_vector(n);

  return A;
}

struct point **allocate_data_points(int n) {
  struct point **data_points = malloc(n * sizeof(struct point *));

  for (int i = 0; i < n; i++)
    data_points[i] = malloc(sizeof(struct point));

  return data_points;
}

void free_matrix(int n, double **A) {
  for (int i = 0; i < n; i++)
    free(A[i]);

  free(A);
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

double *read_vector(int n) {
  double *vector;
  int index;

  vector = allocate_vector(n);

  for (int i = 0; i < n; i++)
    scanf("%d %lf", &index, &vector[i]);

  return vector;
}

double **read_matrix(int n) {
  double **matrix;
  int row_index, column_index;

  matrix = allocate_matrix(n);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      scanf("%d %d %lf", &row_index, &column_index, &matrix[i][j]);

  return matrix;
}

struct point **read_lsp_input(int *n, int *m) {
  struct point **data_points;

  printf("Enter the degree of the polynomial:\n");
  *m = read_size();
  printf("Enter the quantity of data points:\n");
  *n = read_size();
  data_points = allocate_data_points(*n);
  printf("Enter the data points coordinates: t_i y_i\n");

  for (int i = 0; i < *n; i++)
    scanf("%lf %lf", &data_points[i]->t, &data_points[i]->y);

  return data_points;
}

void print_vector(int n, double *x) {
  for (int i = 0; i < n; i++)
    printf("%lf ", x[i]);

  printf("\n");
}

void print_matrix(int n, double **A) {
  for (int i = 0; i < n; i++)
    print_vector(n, A[i]);
}

void print_data_points(int n, struct point **data_points) {
  for (int i = 0; i < n; i++)
    printf("t_%d y_%d: %lf %lf\n", i, i, data_points[i]->t, data_points[i]->y);
}

