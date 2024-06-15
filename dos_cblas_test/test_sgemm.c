#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "cblas.h"

#define max(a,b) \
        ({ __typeof__ (a) _a = (a); \
        __typeof__ (b) _b = (b); \
        _a > _b ? _a : _b; })

// Remember, everything here is COL-MAJOR.
void naive_sgemm_impl(int m, int n, int k,
                float alpha, const float *A, int incRowA, int incColA,
                const float *B, int incRowB, int incColB,
                float beta, float *C, int incRowC, int incColC)
{
    int i, j, d;

    if (beta == 0) {
        if (incRowC == 1) {
            for (i=0; i<n; i++)
                memset(&C[i*incColC], 0, sizeof(float)*m);
        } else {
            for (i=0; i<n; i++)
                for (j=0; j<m; j++)
                    C[i*incColC+j*incRowC] = 0;
        }
    }
    else if (beta != 1) {
        for (i=0; i<n; i++)
            for (j=0; j<m; j++)
                C[i*incColC+j*incRowC] *= beta;
    }

    for (i=0; i<n; i++)
        for (j=0; j<m; j++) 
            for (d=0; d<k; d++)
                C[i*incColC+j*incRowC] += alpha*A[d*incColA+j*incRowA]*B[i*incColB+d*incRowB];
}

void naive_sgemm_(enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB,
                int m, int n, int k,
                float alpha, const float *A, int ldA,
                const float *B, int ldB,
                float beta, float *C, int ldC)
{
//
//  Start the operations.
//
    if (transB==CblasNoTrans) {
        if (transA==CblasNoTrans) {
//
//          Form  C := alpha*A*B + beta*C.
//
            naive_sgemm_impl(m, n, k,
                          alpha,
                          A, 1, ldA,
                          B, 1, ldB,
                          beta,
                          C, 1, ldC);
        } else {
//
//          Form  C := alpha*A**T*B + beta*C
//
            naive_sgemm_impl(m, n, k,
                          alpha,
                          A, ldA, 1,
                          B, 1, ldB,
                          beta,
                          C, 1, ldC);
        }
    } else {
        if (transA==CblasNoTrans) {
//
//          Form  C := alpha*A*B**T + beta*C
//
            naive_sgemm_impl(m, n, k,
                          alpha,
                          A, 1, ldA,
                          B, ldB, 1,
                          beta,
                          C, 1, ldC);
        } else {
//
//          Form  C := alpha*A**T*B**T + beta*C
//
            naive_sgemm_impl(m, n, k,
                          alpha,
                          A, ldA, 1,
                          B, ldB, 1,
                          beta,
                          C, 1, ldC);
        }
    }
}

void naive_sgemm(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE transA, enum CBLAS_TRANSPOSE transB,
                int m, int n, int k,
                float alpha, const float *A, int ldA,
                const float *B, int ldB,
                float beta, float *C, int ldC)
{
    { // Error handling
        int numRowsA;
        int numRowsB;

        if (order==CblasColMajor) {
            numRowsA = (transA==CblasNoTrans) ? m : k;
            numRowsB = (transB==CblasNoTrans) ? k : n;
        } else {
            numRowsB = (transB==CblasNoTrans) ? n : k;
            numRowsA = (transA==CblasNoTrans) ? k : m;
        }

        int info = 0;
        if (order!=CblasColMajor && order!=CblasRowMajor) {
            info = 1;
        } else if (transA!=CblasNoTrans && transA!=CblasTrans
        && transA!=AtlasConj && transA!=CblasConjTrans)
        {
            info = 2;
        } else if (transB!=CblasNoTrans && transB!=CblasTrans
                && transB!=AtlasConj && transB!=CblasConjTrans)
        {
            info = 3;
        } else {
            if (m<0) {
                info = 4;
            } else if (n<0) {
                info = 5;
            } else if (k<0) {
                info = 6;
            } else if (ldA<max(1,numRowsA)) {
                info = 9;
            } else if (ldB<max(1,numRowsB)) {
                info = 11;
            } else if (ldC<max(1,m)) {
                info = 14;
            }
        }
        if (info!=0) {
            printf ("Error! Info: %d\n", info);
            return;
        }
    }

    if (order==CblasColMajor) {
        naive_sgemm_(transA, transB, m, n, k, alpha,
                       A, ldA,
                       B, ldB,
                       beta,
                       C, ldC);
    } else {
        naive_sgemm_(transB, transA, n, m, k, alpha,
                       B, ldB,
                       A, ldA,
                       beta,
                       C, ldC);
    }
}

// ... and everything below is ROW-MAJOR

float* malloc_matrix(uint32_t x,  uint32_t y)
{
    float* m = malloc(x*y*sizeof(float));
    return m;
}

void init_matrix(float* x, uint32_t m, uint32_t n)
{
    uint32_t i;
    for (i = 0; i < m*n; i++)
        x[i] = (float)rand() / RAND_MAX;
}

void init_bias_as_matrix(float* x, uint32_t n, uint32_t d)
{
    uint32_t i;
    for (i = 0; i < d; i++)
        x[i] = (float)rand() / RAND_MAX;
    for (i = 1; i < n; i++)
        memcpy(&x[i*d], x, sizeof(float)*d);
}

double get_dur(clock_t st, clock_t ed)
{
    return 1000.0 * (ed - st) / CLOCKS_PER_SEC;
}

void print_matrix(float* x, int m, int n, const char* name)
{
    printf ("%s:\n", name);
    int i,j;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++)
            printf("%6.2f ", x[i*n+j]);
        printf ("\n");
    }
    printf ("\n");
}

int main(int argc, char** argv)
{
    srand(time(NULL));
    if (argc != 4)
    {
        printf ("Usage:\n");
        printf ("   %s: <n> <in_dim> <out_dim>\n", argv[0]);
        return 0;
    }
    uint32_t n = atoi(argv[1]);
    uint32_t id = atoi(argv[2]);
    uint32_t od = atoi(argv[3]);

    printf ("LinearWithBias test:\n");
    printf ("Input:  [%d x %d]\n", n, id);
    printf ("Weight: [%d x %d]\n", id, od);
    printf ("Output: [%d x %d]\n", n, od);
    printf ("\n");

    clock_t st, ed;
    float* x = malloc_matrix(n, id);    // input
    float* w = malloc_matrix(id, od);   // weight
    float* y = malloc_matrix(n, od);    // output
    ed = clock();
    printf ("Memory alloc: %.2f\n", get_dur(st, ed));

    st = clock();
    init_matrix(x, n, id);
    ed = clock();
    printf ("Init Input: %.2f\n", get_dur(st, ed));

    st = clock();
    init_matrix(w, id, od);
    ed = clock();
    printf ("Init Weight: %.2f\n", get_dur(st, ed));

    st = clock();
    init_bias_as_matrix(y, n, od);
    ed = clock();
    printf ("Init Bias (as output matrix): %.2f\n", get_dur(st, ed));
    printf ("\n");

    if (n <= 5 && id <= 5 && od <= 5) {
        print_matrix(x, n, id, "Input");
        print_matrix(w, id, od, "Weight");
        print_matrix(y, 1, od, "Bias");
    }

    // This is for second test
    float* z = malloc_matrix(n, od);
    memcpy(z, y, sizeof(float)*n*od);
    
    st = clock();
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, od, id, 1.0, x, id, w, od, 1.0, y, od);
    ed = clock();
    printf ("ulmBLAS sgemm: %.2f\n", get_dur(st, ed));

    st = clock();
    naive_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, od, id, 1.0, x, id, w, od, 1.0, z, od);
    ed = clock();
    printf ("naive sgemm: %.2f\n", get_dur(st, ed));

    if (n <= 5 && id <= 5 && od <= 5) {
        print_matrix(y, n, od, "ulmBLAS result");
        print_matrix(z, n, od, "naive result");
    }

    int i,j;
    float delta = 0;
    for (i=0; i<n; i++)
        for (j=0; j<od; j++)
            delta += fabsf(y[i*od+j] - z[i*od+j]);
    printf ("Result delta: %.2f\n", delta/(n*od));

    free(x);
    free(w);
    free(y);
    free(z);
    return 0;
}