#ifndef CFSS_H_
#define CFSS_H_

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
#include <immintrin.h>

#define N 100
#define K 4
#define SEED 2
#define E (K * N - (K * (K + 1)) / 2)
#define R (1 + (((E > N ? E : N) - 1) / 128))

#define FORWARD 40
#define DAYAHEAD 80
#define GAMMA 1.8
//#define KAPPA(S) (pow(S, GAMMA))
#include "k.i"
#define KAPPA(S) (k[S])
#define X(v, i) ((v)[2 * (i)])
#define Y(v, i) ((v)[2 * (i) + 1])

typedef uint32_t agent;
typedef uint32_t edge;
typedef float value;

#include "random.h"
#include "read.h"

typedef struct {
	value p[N * TS], t[N], m[N];
	agent a[2 * (E + 1)], n[2 * N + 1];
	agent s[2 * N], cs[N];
	edge g[N * N];
} stack;

#endif /* CFSS_H_ */
