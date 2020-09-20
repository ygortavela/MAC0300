#include "ep1.h"

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

void matrix_vector_by_row(int n, double **A, double *x, double *b) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      b[i] = b[i] + A[i][j]*x[j];
  }
}

void matrix_vector_by_column(int n, double **A, double *x, double *b) {
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++)
      b[i] = b[i] + A[i][j]*x[j];
  }
}

void matrix_matrix_jki(int n, double **A, double **X, double **B) {
  for (int j = 0; j < n; j++)
    for (int k = 0; k < n; k++)
      for (int i = 0; i < n; i++)
        B[i][j] = B[i][j] + A[i][k]*X[k][j];

}
void matrix_matrix_ikj(int n, double **A, double **X, double **B) {
  for (int i = 0; i < n; i++)
    for (int k = 0; k < n; k++)
      for (int j = 0; j < n; j++)
        B[i][j] = B[i][j] + A[i][k]*X[k][j];
}

void case_one() {
  int n;
  double *x, *y, result;

  n = read_size();
  x = read_vector(n);

  n = read_size();
  y = read_vector(n);


  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  result = dot_product(n, x, y);
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (DEBUG) printf("%lf\n", result);

  free(x);
  free(y);
}

void case_two() {
  int n;
  double *x, result;

  n = read_size();
  x = read_vector(n);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  result = euclidean_norm(n, x);
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (DEBUG) printf("%lf\n", result);

  free(x);
}

void case_three() {
  int n;
  double **A, *x, *b;

  n = read_size();
  A = read_matrix(n);

  n = read_size();
  x = read_vector(n);

  b = allocate_vector(n);
  initialize_vector(n, b);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  matrix_vector_by_row(n, A, x, b);
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (DEBUG) print_vector(n, b);

  free_matrix(n, A);
  free(x);
  free(b);
}

void case_four() {
  int n;
  double **A, *x, *b;

  n = read_size();
  A = read_matrix(n);

  n = read_size();
  x = read_vector(n);

  b = allocate_vector(n);
  initialize_vector(n, b);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  matrix_vector_by_column(n, A, x, b);
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (DEBUG) print_vector(n, b);

  free_matrix(n, A);
  free(x);
  free(b);
}

void case_five() {
  int n;
  double **A, **X, **B;

  n = read_size();
  A = read_matrix(n);

  n = read_size();
  X = read_matrix(n);

  B = allocate_matrix(n);
  initialize_matrix(n, B);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  matrix_matrix_jki(n, A, X, B);
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (DEBUG) print_matrix(n, B);

  free_matrix(n, A);
  free_matrix(n, B);
  free_matrix(n, X);
}

void case_six() {
  int n;
  double **A, **X, **B;

  n = read_size();
  A = read_matrix(n);

  n = read_size();
  X = read_matrix(n);

  B = allocate_matrix(n);
  initialize_matrix(n, B);

  clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
  matrix_matrix_ikj(n, A, X, B);
  clock_gettime(CLOCK_MONOTONIC, &timer.t_end);

  if (DEBUG) print_matrix(n, B);

  free_matrix(n, A);
  free_matrix(n, B);
  free_matrix(n, X);
}

int main(int argc, char* argv[]) {
  int option = argc == 2 ? atoi(argv[1]) : -1;

  if (argc != 2 || option <= 0 || option > 6) {
    printf("Execute ./ep1 #numero_operacao\n");
    printf("Sendo #numero_operacao dado por:\n");
    printf("1 - produto interno entre dois vetores de entrada x e y com tamanho n\n");
    printf("2 - norma euclidiana de um vetor x de tamanho n\n");
    printf("3 - produto matriz A por vetor x seguindo loop da matriz por linhas\n");
    printf("4 - produto matriz A por vetor x seguindo loop da matriz por colunas\n");
    printf("5 - produto matriz A por matriz X seguindo ordem de loop jki\n");
    printf("6 - produto matriz A por matriz X seguindo ordem de loop ikj\n");
  }

  switch (option) {
    case 1:
      case_one();
      break;
    case 2:
      case_two();
      break;
    case 3:
      case_three();
      break;
    case 4:
      case_four();
      break;
    case 5:
      case_five();
      break;
    case 6:
      case_six();
      break;
  }

  if (TIME_BENCHMARK) {
    printf("elapsed time to execute operation: %f\n",
      (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
      (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
  }
}

