#include "util.h"

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

void free_matrix(int n, double **matrix) {
  for (int i = 0; i < n; i++)
    free(matrix[i]);

  free(matrix);
}

int read_size() {
  int n;

  scanf("%d", &n);

  return n;
}

double *read_vector(int n) {
  double *vector;

  vector = allocate_vector(n);

  for (int i = 0; i < n; i++)
    scanf("%lf", &vector[i]);

  return vector;
}

double **read_matrix(int n) {
  double **matrix;

  matrix = allocate_matrix(n);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      scanf("%lf", &matrix[i][j]);

  return matrix;
}

void print_vector(int n, double *x) {
  for (int i = 0; i < n; i++)
    printf("%lf ", x[i]);

  printf("\n");
}

void print_matrix(int n, double **matrix) {
  for (int i = 0; i < n; i++)
    print_vector(n, matrix[i]);
}

double dot_product(int n, double *x, double *y) {
  double dot_product = 0;

  for (int i = 0; i < n; i++)
    dot_product += x[i]*y[i];

  return dot_product;
}

double euclidean_norm(int n, double *x) {
  double normalized_component, normalized_dot_sum = 0, largest_component = largest_vector_component(n, x);

  for (int i = 0; i < n; i++) {
    normalized_component = x[i]/largest_component;
    normalized_dot_sum += normalized_component * normalized_component;
  }

  return largest_component * sqrt(normalized_dot_sum);
}

void matrix_vector_product(int n, double **A, double *x, double *b) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      b[i] = b[i] + A[i][j]*x[j];
  }
}

int cholrow(int n, double **A) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      for (int k = 0; k < j; k++)
        A[i][j] -= A[i][k] * A[j][k];

      A[i][j] /= A[j][j];
    }

    for (int k = 0; k < i; k++)
      A[i][i] -= A[i][k] * A[i][k];

    if (A[i][i] <= 0) return -1;

    A[i][i] = sqrt(A[i][i]);
  }

  return 0;
}

int forwrow(int n, double **A, double *b) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++)
      b[i] -= A[i][j] * b[j];

    if (fabs(A[i][i]) < EPSILON) return -1;

    b[i] /= A[i][i];
  }

  return 0;
}

int backrow(int n, double **A, double *b) {
  for (int j = n - 1; j >= 0; j--) {
    if (fabs(A[j][j]) < EPSILON) return -1;

    b[j] /= A[j][j];

    for (int i = j - 1; i >= 0; i--)
      b[i] -= A[j][i] * b[j];
  }

  return 0;
}

double *deep_copy_vector(int n, double *target_vector) {
  double *vector_copy = malloc(n * sizeof(double));

  for (int i = 0; i < n; i++)
    vector_copy[i] = target_vector[i];

  return vector_copy;
}

void scaled_vector_sum(int n, double lambda, double *u, double *v, double *result) {
  for (int i = 0; i < n; i++)
    result[i] = u[i] + (lambda * v[i]);
}

void scaled_vector_subtraction(int n, double lambda, double *u, double *v, double *result) {
  for (int i = 0; i < n; i++)
    result[i] = u[i] - (lambda * v[i]);
}
