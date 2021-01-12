#include "ep4.h"

void cholesky_method() {
  int n, error;
  double *b, **A;

  n = read_size();
  A = read_matrix(n);
  b = read_vector(n);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  error = cholrow(n, A);
  error = !error ? forwrow(n, A, b) : -1;
  error = !error ? backrow(n, A, b) : -1;
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (!error && TIME_BENCHMARK) {
    printf("Elapsed time to solve Ax = b: %lf\n",
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

void conjugate_gradient_iteration() {
  int n, error;
  double *b, **A;

  n = read_size();
  A = read_matrix(n);
  b = read_vector(n);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  error = cholrow(n, A);
  error = !error ? forwrow(n, A, b) : -1;
  error = !error ? backrow(n, A, b) : -1;
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (!error && TIME_BENCHMARK) {
    printf("Elapsed time to solve Ax = b: %lf\n",
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
    printf("Execute ./ep4 #operation_number\n");
    printf("Linear system solver for n-dimensional positive definite systems\n");
    printf("1 - Conjugate Gradient Iteration\n");
    printf("2 - Cholesky Decomposition\n");
  }

  switch (option) {
    case 1:
      conjugate_gradient_iteration();
      break;
    case 2:
      cholesky_method();
      break;
  }
}
