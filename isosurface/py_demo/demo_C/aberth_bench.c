/* Compiled Ehrlich-Aberth on a Legendre series: the honest per-fiber cost.
   Mirrors demo_C/route_c.py exactly: deflate, Chebyshev+offset init, fixed iters. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define NITER 12
#define MAXN  16

static void aberth(const double *a, int n, double complex *z) {
    for (int i = 0; i < n; ++i) {
        double t = M_PI * (i + 0.5) / n;
        z[i] = cos(t) + 0.1 * I * sin(t);
    }
    for (int it = 0; it < NITER; ++it) {
        double complex w[MAXN];
        for (int i = 0; i < n; ++i) {
            double complex x = z[i];
            double complex P0 = 1.0, P1 = x, dP0 = 0.0, dP1 = 1.0;
            double complex v = a[0] * P0 + a[1] * P1;
            double complex d = a[1] * dP1;
            for (int k = 1; k < n; ++k) {
                double complex P2  = ((2*k+1) * x * P1 - k * P0) / (k + 1);
                double complex dP2 = ((2*k+1) * (P1 + x * dP1) - k * dP0) / (k + 1);
                v += a[k+1] * P2;  d += a[k+1] * dP2;
                P0 = P1; P1 = P2;  dP0 = dP1; dP1 = dP2;
            }
            w[i] = v / d;
        }
        for (int i = 0; i < n; ++i) {
            double complex S = 0.0;
            for (int j = 0; j < n; ++j) if (j != i) S += 1.0 / (z[i] - z[j]);
            z[i] -= w[i] / (1.0 - w[i] * S);
        }
    }
}

int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    long reps = atol(argv[2]);
    double *A = malloc(sizeof(double) * reps * (n + 1));
    srand(42);
    for (long r = 0; r < reps; ++r)
        for (int k = 0; k <= n; ++k)
            A[r*(n+1)+k] = ((double)rand()/RAND_MAX - 0.5) * pow(0.4, k);
    double complex z[MAXN];
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    double sink = 0.0;
    for (long r = 0; r < reps; ++r) { aberth(&A[r*(n+1)], n, z); sink += creal(z[0]); }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double us = ((t1.tv_sec-t0.tv_sec)*1e9 + (t1.tv_nsec-t0.tv_nsec)) / reps / 1e3;
    double flops = 2.0 * NITER * n * n;
    printf("  n=%2d  %8.4f us/fiber  %10.0f fibers/s  model %6.0f flops -> %7.0f Mflop/s  (sink %.3e)\n",
           n, us, 1e6/us, flops, flops/us, sink);
    free(A);
    return 0;
}
