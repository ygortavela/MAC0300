#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "util.h"

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

