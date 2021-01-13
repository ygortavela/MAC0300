#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef _UTIL_H
#define _UTIL_H

#define EPSILON 1e-6
#define MAX_ITERATIONS 100

void initialize_vector(int n, double *b);

void initialize_matrix(int n, double **A);

double *allocate_vector(int n);

double **allocate_matrix(int n);

void free_matrix(int n, double **matrix);

int read_size();

double *read_vector();

double **read_matrix();

void print_vector(int n, double *x);

void print_matrix(int n, double **matrix);

double dot_product(int n, double *x, double *y);

void matrix_vector_product(int n, double **A, double *x, double *b);

int cholrow(int n, double **A);

int forwrow(int n, double **A, double *b);

int backrow(int n, double **A, double *b);

double *deep_copy_vector(int n, double *target_vector);

void scaled_vector_sum(int n, double lambda, double *u, double *v, double *result);

void scaled_vector_subtraction(int n, double lambda, double *u, double *v, double *result);

#endif
