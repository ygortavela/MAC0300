#include "ep3.h"

void compute_reflector(int n, int k, double *x, double *gamma, double *tau) {
  int max_index = largest_vector_component_index(n, k, x);
  double beta = x[max_index];

  if (fabs(beta) < EPSILON) *gamma = 0;
  else {
    for (int i = k; i < n; i++) x[i] /= beta;

    *tau = euclidean_norm(n, k, x);

    if (x[k] < 0) *tau *= -1;

    x[k] += *tau;
    *gamma = x[k] / (*tau);

    for (int i = k + 1; i < n; i++)
      x[i] /= x[k];

    x[k] = 1;
    *tau *= beta;
  }
}

void compute_QB(int n, int m, int k, double gamma, double **B) {
  double *v_t = allocate_vector(n - k), *temp = allocate_vector(m - k - 1);

  initialize_vector(m - k - 1, temp);

  for (int i = k; i < n; i++)
    v_t[i - k] = gamma * B[k][i];

  for (int j = k + 1; j < m; j++)
    for (int i = k; i < n; i++)
        temp[j - k - 1] += v_t[i - k] * B[j][i];

  for (int j = k + 1; j < m; j++)
    for (int i = k; i < n; i++)
      B[j][i] -= B[k][i] * temp[j - k - 1];

  free(v_t);
  free(temp);
}

int decompose_transpose_to_QR(int n, int m, double **A, double *gamma) {
  double tau;

  for (int k = 0; k < m; k++) {
    compute_reflector(n, k, A[k], &gamma[k], &tau);

    if (gamma[k] == 0) return -1;

    compute_QB(n, m, k, gamma[k], A);
    A[k][k] = - tau;
  }

  return 0;
}

double *apply_reflectors_transpose(int n, int m, struct point **data_points, double **A, double *gamma) {
  double temp, *v_t = allocate_vector(n), *c = allocate_vector(n);

  for (int i = 0; i < n; i++) c[i] = data_points[i]->y;

  for (int k = 0; k < m; k++) {
    temp = 0;
    v_t[k] = gamma[k];

    for (int i = k + 1; i < n; i++)
      v_t[i] = gamma[k] * A[k][i];

    for (int i = k; i < n; i++)
        temp += v_t[i] * c[i];

    c[k] -= temp;

    for (int i = k + 1; i < n; i++)
        c[i] -= A[k][i] * temp;
  }

  free(v_t);

  return c;
}

void full_rank() {
  int n, m, error;
  struct point **data_points;
  double *gamma, **A, *c;

  data_points = read_lsp_input(&n, &m);

  if (m > n) {
    printf("The system is not overdetermined because the degree of the polynomial is greater than the number of given data points.\n");
    free_data_points(n, data_points);
    return;
  }

  A = build_transpose_coefficient_matrix_on_std_basis(n, m, data_points);
  gamma = allocate_vector(m);
  error = decompose_transpose_to_QR(n, m, A, gamma);

  if (!error) {
    c = apply_reflectors_transpose(n, m, data_points, A, gamma);
    backrow(m, A, c);
    print_polynomial_solution(m, c);
    free(c);
  } else {
    printf("Overdetermined system hasn't full rank.");
  }

  free_matrix(m, A);
  free_data_points(n, data_points);
  free(gamma);
}

void compute_reflector_without_scaling(int n, int k, double *x, double *gamma, double *tau) {
  *tau = euclidean_norm(n, k, x);

  if (x[k] < 0) *tau *= -1;

  x[k] += *tau;
  *gamma = x[k] / (*tau);

  for (int i = k + 1; i < n; i++)
    x[i] /= x[k];

  x[k] = 1;
}

int decompose_transpose_to_QR_with_column_interchange(int n, int m, double **A, double *gamma, int *permutation) {
  int k, pivot_index, temp;
  double tau, largest_row, *cached_norms = build_cached_row_norms_vector(m, n, A);
  double matrix_norm = frobenius_norm_using_cached_norms_vector(m, cached_norms);

  pivot_index = largest_vector_component_index(m, 0, cached_norms);
  largest_row = cached_norms[pivot_index];

  for (k = 0; k < m && largest_row > EPSILON * matrix_norm ; k++) {
    if (pivot_index != k) {
      swap_int_vector_value(k, pivot_index, permutation);
      interchange_pivot_row(k, pivot_index, A);
      swap_double_vector_value(k, pivot_index, cached_norms);
    }

    compute_reflector_without_scaling(n, k, A[k], &gamma[k], &tau);
    compute_QB(n, m, k, gamma[k], A);
    A[k][k] = - tau;

    if (k + 1 < m) {
      update_cached_norms_vector(m, k + 1, cached_norms, A);
      pivot_index = largest_vector_component_index(m, k + 1, cached_norms);
      largest_row = cached_norms[pivot_index];
    }
  }

  free(cached_norms);

  return k;
}

void rank_deficient() {
  int n, m, r;
  struct point **data_points;
  double *gamma, **A, *c;
  int *permutation;

  data_points = read_lsp_input(&n, &m);

  if (m > n) {
    printf("The system is not overdetermined because the degree of the polynomial is greater than the number of given data points.\n");
    free_data_points(n, data_points);
    return;
  }

  permutation = malloc(m * sizeof(int));

  for (int i = 0; i < m; i++)
    permutation[i] = i;

  A = build_transpose_coefficient_matrix_on_std_basis(n, m, data_points);
  system_rescale(m, n, A, data_points);
  gamma = allocate_vector(m);
  r = decompose_transpose_to_QR_with_column_interchange(n, m, A, gamma, permutation);
  c = apply_reflectors_transpose(n, r, data_points, A, gamma);
  backrow(r, A, c);

  assemble_permuted_solution(r, m, c, permutation);
  print_polynomial_solution(m, c);

  free_matrix(m, A);
  free_data_points(n, data_points);
  free(gamma);
  free(permutation);
  free(c);
}

int main(int argc, char* argv[]) {
  int option = argc == 2 ? atoi(argv[1]) : -1;

  if (argc != 2 || option <= 0 || option > 2) {
    printf("Execute ./ep3 #operation_number\n");
    printf("Least Squares Problem Solver\n");
    printf("1 - Full rank case\n");
    printf("2 - Rank deficient case\n");
  }

  switch (option) {
    case 1:
      full_rank();
      break;
    case 2:
      rank_deficient();
      break;
  }
}
