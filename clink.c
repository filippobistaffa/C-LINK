#include "clink.h"

struct timeval t1, t2;
stack sol;
value min;

__attribute__((always_inline)) inline
void contract(edge *g, agent *a, const agent *n, agent v1, agent v2, float *link) {

	register agent i, e, f, m = n[N];
	register const agent *p = n + N + 1;

	do if ((i = *(p++)) != v1)
		if ((e = g[i * N + v2])) {
			if ((f = g[i * N + v1])) link[f] = FLT_MAX;
			g[i * N + v1] = g[v1 * N + i] = e;
			X(a, e) = v1;
			Y(a, e) = i;
		}
	while (--m);
}

__attribute__((always_inline)) inline
void merge(stack *st, agent v1, agent v2) {

	register agent a, b, i, j, min = v1, max = v2, *p = st->n + N + 1;

	if (Y(st->s, max) < Y(st->s, min)) {
		b = max;
		max = min;
		min = b;
	}

	a = X(st->s, min);
	b = X(st->s, max);
	max = Y(st->s, max);
	Y(st->s, v1) = min = Y(st->s, min);
	agent c[b];
	X(st->s, v1) = a + b;
	memcpy(c, st->cs + max, sizeof(agent) * b);
	memmove(st->cs + min + a + b, st->cs + min + a, sizeof(agent) * (max - min - a));
	memmove(st->cs + min, st->cs + min, sizeof(agent) * a);
	memcpy(st->cs + min + a, c, sizeof(agent) * b);

	if ((j = st->n[st->n[N] + N]) != v2) {
		st->n[j] = st->n[v2];
		st->n[st->n[v2]] = j;
		st->n[v2] = st->n[N] + N;
	}

	j = --st->n[N];

	do if ((i = *(p++)) != v1) {
		a = Y(st->s, i);
		if (a > min && a < max) Y(st->s, i) = a + b;
	} while (--j);
}

__attribute__((always_inline))
inline void mergeprof(stack *st, agent v1, agent v2) {

	register uint_fast64_t i;
	register __m128 a, b = _mm_set1_ps(FLT_MAX);
	st->t[v1] += st->t[v2];

	for (i = 0; i < TS; i += 4) {
		a = _mm_add_ps(_mm_load_ps(st->p + v1 * TS + i), _mm_load_ps(st->p + v2 * TS + i));
		b = _mm_min_ps(a, b);
		_mm_store_ps(st->p + v1 * TS + i, a);
	}

	b = _mm_min_ps(b, _mm_shuffle_ps(b, b, _MM_SHUFFLE(2, 1, 0, 3)));
	st->m[v1] = _mm_min_ps(b, _mm_shuffle_ps(b, b, _MM_SHUFFLE(1, 0, 3, 2)))[0];
}

__attribute__((always_inline))
inline value val(stack *st) {

	register value h, j, r = 0;
	register agent i, nn = st->n[N];
	register const agent *p = st->n + N + 1;

	do {
		i = *(p++);
		h = st->t[i] - (j = st->m[i] * TS);
		r += h * DAYAHEAD + j * FORWARD + KAPPA(X(st->s, i));
	} while (--nn);

	return r;
}

#define GAIN(MER, CUR) (val(MER) - val(CUR))

void clink(stack *st, float *link, value tot) {

	if (tot < min) { sol = *st; min = tot; }

	float bestlinkage = FLT_MAX;
	agent v1 = 0, v2 = 0;
	edge beste = 0;

	for (edge e = 1; e < E + 1; e++)
		if (link[e] < bestlinkage) {
			bestlinkage = link[e];
			v1 = X(st->a, e);
			v2 = Y(st->a, e);
			beste = e;
		}

	if (bestlinkage < 0) {
		link[beste] = FLT_MAX;
		//printf("will contract %u (%u -- %u) %f\n", beste, v1, v2, bestlinkage);
		merge(st, v1, v2); // merge best coalitions
		mergeprof(st, v1, v2);
		contract(st->g, st->a, st->n, v1, v2, link);
		agent *a = st->n + N + 1;
		agent m = st->n[N];
		edge e = 0;
		do {
			v2 = *(a++);
			if (v1 != v2 && (e = st->g[v1 * N + v2])) {
				register stack mst = *st;
				merge(&mst, v1, v2);
				mergeprof(&mst, v1, v2);
				link[e] = GAIN(&mst, st);
			} else link[e] = FLT_MAX;
		} while (--m);
		clink(st, link, tot + bestlinkage);
	}
}

inline void createedge(edge *g, agent *a, agent v1, agent v2, edge e) {

	g[v1 * N + v2] = g[v2 * N + v1] = e;
	a[e * 2] = v1;
	a[e * 2 + 1] = v2;
}

void createScaleFree(edge *g, agent *a) {

	uint_fast8_t deg[N] = {0};
	register uint_fast64_t d, i, j, h, k = 1, q, t = 0;
	register int p;

	for (i = 1; i <= K; i++) {
		for (j = 0; j < i; j++) {
			createedge(g, a, i, j, k++);
			deg[i]++;
			deg[j]++;
		}
	}

	for (i = K + 1; i < N; i++) {
		t &= ~((1UL << i) - 1);
		for (j = 0; j < K; j++) {
			d = 0;
			for (h = 0; h < i; h++)
				if (!((t >> h) & 1)) d += deg[h];
			if (d > 0) {
				p = nextInt(d);
				q = 0;
				while (p >= 0) {
					if (!((t >> q) & 1)) p = p - deg[q];
					q++;
				}
				q--;
				t |= 1UL << q;
				createedge(g, a, i, q, k++);
				deg[i]++;
				deg[q]++;
			}
		}
	}
}

void printcs(stack *st) {

	register const agent *p = st->n + N + 1;
        register agent i, j, m = st->n[N];

	do {
		i = *(p++);
                printf("{ ");
                for (j = 0; j < X(st->s, i); j++)
                	printf("%s%u%s ", i == st->cs[Y(st->s, i) + j] ? "<" : "", 
			       st->cs[Y(st->s, i) + j], i == st->cs[Y(st->s, i) + j] ? ">" : "");
                printf("}\n");
        } while (--m);
}

int main(int argc, char *argv[]) {

	stack *st = (stack *)malloc(sizeof(stack));
	memset(st->g, 0, sizeof(edge) * N * N);
	st->n[N] = N;

	for (int i = 0; i < N; i++) {
		X(st->s, i) = 1;
		Y(st->s, i) = st->cs[i] = i;
		st->n[st->n[i] = N + i + 1] = i;
	}

	init(SEED);
	createScaleFree(st->g, st->a);
	init(SEED);
	read(st->p, st->t, st->m);

	gettimeofday(&t1, NULL);
	float *link = (float *)malloc(sizeof(float) * (E + 1));
	//#pragma omp parallel for schedule(dynamic)
	for (edge e = 1; e < E + 1; e++) {
		agent v1 = X(st->a, e);
		agent v2 = Y(st->a, e);
		register stack mst = *st;
		merge(&mst, v1, v2);
		mergeprof(&mst, v1, v2);
		link[e] = GAIN(&mst, st);
	}

	sol = *st;
	value in = min = val(st);
	clink(st, link, in);
	free(st);
	gettimeofday(&t2, NULL);
	//printcs(&sol);
	printf("%u,%u,%f,%f,%f,%f\n", N, SEED, in, min, (in - min) / in,
	       (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);
	return 0;
}
