#include "ep4.h"

void conjugate_gradient_iteration() {
  int n;
  double *b, **A, *x, *previous_residual, *residual, *search_direction, *streched_search_direction;
  double alfa, beta, residual_dot_product;

  read_system_data(&n, A, b);
  x = allocate_vector(n);
  initialize_vector(n, x);
  residual = allocate_vector(n);
  previous_residual = deep_copy_vector(n, b);
  search_direction = deep_copy_vector(n, b);
  streched_search_direction = allocate_vector(n);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);

  for (int i = 0; i < MAX_ITERATIONS; i++) {
    matrix_vector_product(n, A, search_direction, streched_search_direction);
    alfa = dot_product(n, previous_residual, previous_residual)/dot_product(n, search_direction, streched_search_direction);
    scaled_vector_sum(n, alfa, x, previous_residual, x);
    scaled_vector_subtraction(n, alfa, previous_residual, streched_search_direction, residual);
    residual_dot_product = dot_product(n, residual, residual);

    if (fabs(residual_dot_product) > EPSILON) break;

    beta = residual_dot_product/dot_product(n, previous_residual, previous_residual);
    scaled_vector_sum(n, beta, residual, search_direction, search_direction);
  }

  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (TIME_BENCHMARK) {
    printf("Elapsed time to solve Ax = b: %lf\n",
      (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
      (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
  }

  if (DEBUG)
    print_vector(n, b);

  free_matrix(n, A);
  free(b);
  free(x);
  free(previous_residual);
  free(residual);
  free(search_direction);
  free(streched_search_direction);
}

void cholesky_method() {
  int n, error;
  double *b, **A;

  read_system_data(&n, A, b);

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
