#ifndef _UTIL_H
#define _UTIL_H

int pivot_row_index(int n, double **A, int from_index);

void interchange_pivot_row(int k, int pivot_index, double **A);

void initialize_vector(int n, double *b);

void initialize_matrix(int n, double **A);

double *allocate_vector(int n);

double **allocate_matrix(int n);

void free_matrix(int n, double **A);

int read_size();

double *read_vector();

double **read_matrix();

void print_vector();

void print_matrix();

#endif
