#include "ep2.h"

int cholcol(int n, double **A) {
  for (int j = 0; j < n; j++) {
    if (A[j][j] <= 0) return -1;

    A[j][j] = sqrt(A[j][j]);

    for (int i = j + 1; i < n; i++)
      A[i][j] /= A[j][j];

    for (int k = j + 1; k < n; k++)
      for (int i = k; i < n; i++)
        A[i][k] -= A[i][j] * A[k][j];
  }

  return 0;
}

int forwcol(int n, double **A, double *b) {
  for (int j = 0; j < n; j++) {
    if (fabs(A[j][j]) < EPSILON) return -1;

    b[j] /= A[j][j];

    for (int i = j + 1; i < n; i++)
      b[i] -= A[i][j] * b[j];
  }

  return 0;
}

int backcol(int n, double **A, double *b, int trans) {
  if (trans == 0) {
    for (int j = n - 1; j >= 0; j--) {
      if (fabs(A[j][j]) < EPSILON) return -1;

      b[j] /= A[j][j];

      for (int i = j - 1; i >= 0; i--)
        b[i] -= A[i][j] * b[j];
    }
  } else if (trans == 1) {
    for (int i = n - 1; i >= 0; i--) {
      for (int j = n - 1; j > i; j--)
        b[i] -= A[j][i] * b[j];

      if (fabs(A[i][i]) < EPSILON) return -1;

      b[i] /= A[i][i];
    }
  }

  return 0;
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

int backrow(int n, double **A, double *b, int trans) {
  if (trans == 0) {
    for (int i = n - 1; i >= 0; i--) {
      for (int j = n - 1; j > i; j--)
        b[i] -= A[i][j] * b[j];

      if (fabs(A[i][i]) < EPSILON) return -1;

      b[i] /= A[i][i];
    }
  } else if (trans == 1) {
    for (int j = n - 1; j >= 0; j--) {
      if (fabs(A[j][j]) < EPSILON) return -1;

      b[j] /= A[j][j];

      for (int i = j - 1; i >= 0; i--)
        b[i] -= A[j][i] * b[j];
    }
  }

  return 0;
}

void cholesky_method(int option) {
  int n, error;
  double *b, **A;

  n = read_size();
  A = read_matrix(n);
  b = read_vector(n);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  if (option == COLUMN) error = cholcol(n, A);
  else if (option == ROW) error = cholrow(n, A);
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (!error && TIME_BENCHMARK) {
    printf("Elapsed time to solve A = GG': %lf\n",
      (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
      (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
  }

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  if (option == COLUMN) error = !error ? forwcol(n, A, b) : -1;
  else if (option == ROW) error = !error ? forwrow(n, A, b) : -1;
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (!error && TIME_BENCHMARK) {
    printf("Elapsed time to solve Gy = b': %lf\n",
      (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
      (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
  }

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  if (option == COLUMN) error = !error ? backcol(n, A, b, 1) : -1;
  else if (option == ROW) error = !error ? backrow(n, A, b, 1) : -1;
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (!error && TIME_BENCHMARK) {
    printf("Elapsed time to solve G'x = y: %lf\n",
      (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
      (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
  }

  if (DEBUG) {
    if (!error) print_vector(n, b);
    else printf("Linear system has infinite solutions\n");
  }

  free_matrix(n, A);
  free(b);
}

int main(int argc, char* argv[]) {
  int option = argc == 2 ? atoi(argv[1]) : -1;

  if (argc != 2 || option <= 0 || option > 2) {
    printf("Execute ./ep2 #operation_number\n");
    printf("Linear system solver of lenght n\n");
    printf("1 - Positive definite system using Cholesky Decomposition method by COLUMNS\n");
    printf("2 - Positive definite system using Cholesky Decomposition method by ROWS\n");
  }

  switch (option) {
    case 1:
      cholesky_method(COLUMN);
      break;
    case 2:
      cholesky_method(ROW);
      break;
  }
}

