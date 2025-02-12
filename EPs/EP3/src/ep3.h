#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "util.h"

#ifndef _EP3_H
#define _EP3_H

void compute_reflector(int n, int k, double *x, double *gamma, double *tau);

void compute_QB(int n, int m, int k, double gamma, double **B);

int decompose_transpose_to_QR(int n, int m, double **A, double *gamma);

double *apply_reflectors_transpose(int n, int m, struct point **data_points, double **A, double *gamma);

void full_rank();

void compute_reflector_without_scaling(int n, int k, double *x, double *gamma, double *tau);

int decompose_transpose_to_QR_with_column_interchange(int n, int m, double **A, double *gamma, int *permutation);

void rank_deficient();

#endif
