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

#ifndef _EP3_H
#define _EP3_H

struct timer_info {
  struct timespec t_start;
  struct timespec t_end;
};

struct timer_info timer;

void compute_reflector(int n, int k, double *x, double *gamma, double *tau);

void compute_QB(int n, int m, int k, double gamma, double **B);

int decompose_transpose_to_QR(int n, int m, double **A, double *gamma);

double *apply_reflectors(int n, int m, struct point **data_points, double **A, double *gamma);

void full_rank();


#endif
