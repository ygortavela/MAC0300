#include <time.h>

#include "util.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef TIME_BENCHMARK
#define TIME_BENCHMARK 0
#endif

#ifndef _EP4_H
#define _EP4_H

struct timer_info {
    struct timespec t_start;
    struct timespec t_end;
};

struct timer_info timer;

void conjugate_gradient_iteration();

void cholesky_method();

#endif
