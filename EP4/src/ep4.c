#include "ep4.h"

void conjugate_gradient_iteration() {
  int n, i;
  double *b, **A, *x_previous, *x, *previous_residual, *residual, *previous_search_direction, *search_direction, *streched_search_direction;
  double alfa, beta, previous_residual_dot_product, residual_dot_product;

  n = read_size();
  A = read_matrix(n);
  b = read_vector(n);
  x_previous = allocate_vector(n);
  initialize_vector(n, x_previous);
  previous_residual = deep_copy_vector(n, b);
  previous_search_direction = deep_copy_vector(n, b);
  streched_search_direction = allocate_vector(n);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);

  for (i = 0; i < MAX_ITERATIONS; i++) {
    x = allocate_vector(n);
    residual = allocate_vector(n);
    search_direction = allocate_vector(n);
    initialize_vector(n, streched_search_direction);
    matrix_vector_product(n, A, previous_search_direction, streched_search_direction);
    previous_residual_dot_product = dot_product(n, previous_residual, previous_residual);

    alfa = previous_residual_dot_product/dot_product(n, previous_search_direction, streched_search_direction);
    scaled_vector_sum(n, alfa, x_previous, previous_residual, x);
    scaled_vector_subtraction(n, alfa, previous_residual, streched_search_direction, residual);
    residual_dot_product = dot_product(n, residual, residual);

    if (residual_dot_product < EPSILON) {
      free(x_previous);
      free(previous_residual);
      free(previous_search_direction);
      break;
    }

    beta = residual_dot_product/previous_residual_dot_product;
    scaled_vector_sum(n, beta, residual, previous_search_direction, search_direction);

    free(x_previous);
    free(previous_residual);
    free(previous_search_direction);
    x_previous = x;
    previous_residual = residual;
    previous_search_direction = search_direction;
  }

  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  printf("Number of iterations: %d\n", i);

  if (TIME_BENCHMARK) {
    printf("Elapsed time to solve Ax = b: %lf\n",
      (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
      (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
  }

  if (DEBUG) print_vector(n, x);

  free_matrix(n, A);
  free(b);
  free(x);
  free(residual);
  free(search_direction);
  free(streched_search_direction);
}

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
