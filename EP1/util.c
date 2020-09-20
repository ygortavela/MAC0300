#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"

double largest_vector_component(int n, double *x) {
  double temp, max = fabs(x[0]);

  for (int i = 1; i < n; i++) {
    temp = fabs(x[i]);

    if (temp > max)
      max = temp;
  }

  return max;
}

void initialize_vector(int n, double *b) {
  for (int i = 0; i < n; i++)
    b[i] = 0;
}

void initialize_matrix(int n, double **A) {
  for (int i = 0; i < n; i++)
    initialize_vector(n, A[i]);
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

void free_matrix(int n, double **A) {
  for (int i = 0; i < n; i++)
    free(A[i]);

  free(A);
}

double *read_vector() {
  int n;
  double *vector;

  scanf("%d", &n);
  vector = allocate_vector(n);

  for (int i = 0; i < n; i++)
    scanf("%lf", vector[i]);

  return vector;
}

double **read_matrix() {
  int n;
  double **matrix;

  scanf("%d", &n);
  matrix = allocate_matrix(n);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      scanf("%lf", &matrix[i][j]);

  return matrix;
}
