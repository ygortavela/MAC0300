#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "util.h"

double largest_vector_component(int n, double *x) {
  double temp, max = fabs(x[0]);

  for (int i = 1; i < n; i++) {
    temp = fabs(x[i]);

    if (temp > max)
      max = temp;
  }

  return max;
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

void free_matrix(int n, double **A) {
  for (int i = 0; i < n; i++)
    free(A[i]);

  free(A);
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

void print_vector(int n, double *x) {
  for (int i = 0; i < n; i++)
    printf("%lf ", x[i]);

  printf("\n");
}

void print_matrix(int n, double **A) {
  for (int i = 0; i < n; i++)
    print_vector(n, A[i]);
}

