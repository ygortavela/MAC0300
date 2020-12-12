#include "ep3.h"

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

int lucol(int n, double **A, int *p) {
  int pivot_index;

  for (int k = 0; k < n - 1; k++) {
    pivot_index = pivot_row_index(n, A, k);

    if (fabs(A[pivot_index][k]) < EPSILON) {
      p[k] = -1;

      return -1;
    } else {
      p[k] = pivot_index;

      if (k != pivot_index) interchange_pivot_row(k, pivot_index, A);

      for (int i = k + 1; i < n; i++)
        A[i][k] /= A[k][k];

      for (int j = k + 1; j < n; j++)
        for (int i = k + 1; i < n; i++)
          A[i][j] -= A[i][k] * A[k][j];
    }
  }

  if (A[n - 1][n - 1] == 0) {
    p[n - 1] = -1;

    return -1;
  }

  p[n - 1] = n - 1;

  return 0;
}

int sscol(int n, double **A, int *p, double *b) {
  double temp;

  for (int k = 0; k < n - 1; k++) {
    temp = b[p[k]];
    b[p[k]] = b[k];
    b[k] = temp;
  }

  for (int j = 0; j < n - 1; j++)
    for (int i = j + 1; i < n; i++)
      b[i] -= A[i][j] * b[j];

  return backcol(n, A, b, 0);
}

int lurow(int n, double **A, int *p) {
  int pivot_index;

  for (int k = 0; k < n - 1; k++) {
    pivot_index = pivot_row_index(n, A, k);

    if (fabs(A[pivot_index][k]) < EPSILON) {
      p[k] = -1;

      return -1;
    } else {
      p[k] = pivot_index;

      if (k != pivot_index) interchange_pivot_row(k, pivot_index, A);

      for (int i = k + 1; i < n; i++)
        A[i][k] /= A[k][k];

      for (int i = k + 1; i < n; i++)
        for (int j = k + 1; j < n; j++)
          A[i][j] -= A[i][k] * A[k][j];
    }
  }

  if (A[n - 1][n - 1] == 0) {
    p[n - 1] = -1;

    return -1;
  }

  p[n - 1] = n - 1;

  return 0;
}

int ssrow(int n, double **A, int *p, double *b) {
  double temp;

  for (int k = 0; k < n - 1; k++) {
    temp = b[p[k]];
    b[p[k]] = b[k];
    b[k] = temp;
  }

  for (int i = 1; i < n; i++)
    for (int j = 0; j < i; j++)
      b[i] -= A[i][j] * b[j];

  return backrow(n, A, b, 0);
}

/*
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

void gaussian_method(int option) {
  int n, error;
  double *b, **A;
  int *p;

  n = read_size();
  A = read_matrix(n);
  b = read_vector(n);
  p = malloc(n * sizeof(int));

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  if (option == COLUMN) error = lucol(n, A, p);
  else if (option == ROW) error = lurow(n, A, p);
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (!error && TIME_BENCHMARK) {
    printf("Elapsed time to solve PA = LU: %lf\n",
      (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
      (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
  }

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  if (option == COLUMN) error = !error ? sscol(n, A, p, b) : -1;
  else if (option == ROW) error = !error ? ssrow(n, A, p, b) : -1;
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (!error && TIME_BENCHMARK) {
    printf("Elapsed time to solve LUx = Pb: %lf\n",
      (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
      (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
  }

  if (DEBUG) {
    if (!error) print_vector(n, b);
    else printf("Linear system has infinite solutions\n");
  }

  free_matrix(n, A);
  free(b);
  free(p);
}
*/

double compute_reflector_gamma(int n, int k, double *x) {
  double gamma, tau, beta = largest_vector_component(n, k, x);

  if (fabs(beta) < EPSILON) gamma = 0;
  else {
    for (int i = k; i < n; i++) x[i] /= beta;

    tau = euclidean_norm(n, k, x);

    if (x[k] < 0) tau = -tau;

    x[k] += tau;
    gamma = x[k]/tau;

    for (int i = k + 1; i < n; i++)
      x[i] /= x[k];

    x[k] = tau * beta;
  }

  return gamma;
}

void compute_QB(int n, int m, int k, double gamma, double *u, double **B) {
  double *v_t = allocate_vector(n - k), *temp = allocate_vector(m - k - 1);

  v_t[0] = gamma;

  for (int i = k + 1; i < n; i++)
    v_t[i - k - 1] = gamma * u[i];

  for (int j = k + 1; j < m; j++)
    for (int i = k; i < n; i++)
        temp[j - k - 1] += v_t[i - k] * B[j][i];

  for (int j = k + 1; j < m; j++)
    for (int i = k; i < n; i++)
        B[j][i] -= u[i] * temp[j - k - 1];

  free(v_t);
  free(temp);
}


void full_rank() {
  int n, m;
  struct point **data_points;
  double **A;

  data_points = read_lsp_input(&n, &m);

  A = build_transpose_coefficient_matrix_on_std_basis(n, m, data_points);
  print_matrix(m, n, A);

  free_matrix(m, A);
  free_data_points(n, data_points);
}

int main(int argc, char* argv[]) {
  int option = argc == 2 ? atoi(argv[1]) : -1;

  if (argc != 2 || option <= 0 || option > 4) {
    printf("Execute ./ep3 #operation_number\n");
    printf("Least Squares Problem Solver\n");
    printf("1 - Full rank case (if the program detect that the system is rank deficient, it will try to solve using the rank deficient method)\n");
    printf("2 - Rank deficient case\n");
  }

  switch (option) {
    case 1:
      full_rank();
      break;
    case 2:
      break;
  }
}

