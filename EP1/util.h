#ifndef _UTIL_H
#define _UTIL_H

double largest_vector_component(int n, double *x);

void initialize_vector(int n, double *b);

void initialize_matrix(int n, double **A);

double *allocate_vector(int n);

double **allocate_matrix(int n);

void free_matrix(int n, double **A);

double *read_vector();

double **read_matrix();

#endif
