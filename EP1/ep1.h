#ifndef _EP1_H
#define _EP1_H

double dot_product(int n, double *x, double *y);

double euclidean_norm(int n, double *x);

void matrix_vector_by_row(int n, double **A, double *x, double *b);
void matrix_vector_by_column(int n, double **A, double *x, double *b);

void matrix_matrix_ijk(int n, double **A, double **X, double **B);
void matrix_matrix_ikj(int n, double **A, double **X, double **B);

void case_one();
void case_two();
void case_three();
void case_four();
void case_five();
void case_six();

#endif
