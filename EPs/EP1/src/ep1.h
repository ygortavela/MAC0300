#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "util.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef TIME_BENCHMARK
#define TIME_BENCHMARK 0
#endif

#ifndef _EP1_H
#define _EP1_H
struct timer_info {
    struct timespec t_start;
    struct timespec t_end;
};

struct timer_info timer;

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
