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

#define EPSILON 1e-6
#define COLUMN 1
#define ROW 2

struct timer_info {
    struct timespec t_start;
    struct timespec t_end;
};

struct timer_info timer;

int cholcol(int n, double **A);
int forwcol(int n, double **A, double *b);
int backcol(int n, double **A, double *b, int trans);

int cholrow(int n, double **A);
int forwrow(int n, double **A, double *b);
int backrow(int n, double **A, double *b, int trans);

#endif
