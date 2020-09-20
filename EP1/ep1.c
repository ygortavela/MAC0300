#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "util.h"

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
  initialize_vector(n, b);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      b[i] = b[i] + A[i][j]*x[j];
  }
}

void matrix_vector_by_column(int n, double **A, double *x, double *b) {
  initialize_vector(n, b);

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++)
      b[i] = b[i] + A[i][j]*x[j];
  }
}

void matrix_matrix_ijk(int n, double **A, double **X, double **B) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        B[i][j] = B[i][j] + A[i][k]*X[k][j];
  }
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
    printf("5 - produto matriz A por matriz X seguindo ordem de loop ijk\n");
    printf("6 - produto matriz A por matriz X seguindo ordem de loop ikj\n");
  }

  printf("%d\n", option);
}

