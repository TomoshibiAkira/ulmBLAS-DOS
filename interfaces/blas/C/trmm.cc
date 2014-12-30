#include BLAS_HEADER
#include <algorithm>
#include <cctype>
#include <cmath>
#include <interfaces/blas/C/transpose.h>
#include <interfaces/blas/C/xerbla.h>
#include <ulmblas/level3/trlmm.h>
#include <ulmblas/level3/trumm.h>

#include <ulmblas/auxiliary/printmatrix.h>

extern "C" {

void
ULMBLAS(dtrmm)(enum CBLAS_SIDE      side,
               enum CBLAS_UPLO      upLo,
               enum CBLAS_TRANSPOSE transA,
               enum CBLAS_DIAG      diag,
               int                  m,
               int                  n,
               double               alpha,
               const double         *A,
               int                  ldA,
               double               *B,
               int                  ldB)
{
    bool left     = (side==CblasLeft);
    bool lower    = (upLo==CblasLower);
    bool trans    = (transA==CblasTrans || transA==CblasConjTrans);
    bool unitDiag = (diag==CblasUnit);

//
//  Start the operations.
//
    if (left) {
        if (lower) {
            if (!trans) {
                ulmBLAS::trlmm(m, n, alpha, unitDiag, A, 1, ldA, B, 1, ldB);
            } else {
                ulmBLAS::trumm(m, n, alpha, unitDiag, A, ldA, 1, B, 1, ldB);
            }
        } else {
            if (!trans) {
                ulmBLAS::trumm(m, n, alpha, unitDiag, A, 1, ldA, B, 1, ldB);
            } else {
                ulmBLAS::trlmm(m, n, alpha, unitDiag, A, ldA, 1, B, 1, ldB);
            }
        }
    } else {
        if (lower) {
            if (!trans) {
                ulmBLAS::trumm(n, m, alpha, unitDiag, A, ldA, 1, B, ldB, 1);
            } else {
                ulmBLAS::trlmm(n, m, alpha, unitDiag, A, 1, ldA, B, ldB, 1);
            }
        } else {
            if (!trans) {
                ulmBLAS::trlmm(n, m, alpha, unitDiag, A, ldA, 1, B, ldB, 1);
            } else {
                ulmBLAS::trumm(n, m, alpha, unitDiag, A, 1, ldA, B, ldB, 1);
            }
        }
    }
}

void
CBLAS(dtrmm)(enum CBLAS_ORDER     order,
             enum CBLAS_SIDE      side,
             enum CBLAS_UPLO      upLo,
             enum CBLAS_TRANSPOSE transA,
             enum CBLAS_DIAG      diag,
             int                  m,
             int                  n,
             double               alpha,
             const double         *A,
             int                  ldA,
             double               *B,
             const int            ldB)
{
    bool left     = (side==CblasLeft);

//
//  Set  numRowsA and numRowsB as the number of rows of A and B
//
    int numRowsA = (left) ? m : n;

//
//  Test the input parameters
//
    int info = 0;

    if (order!=CblasColMajor && order!=CblasRowMajor) {
        info = 1;
    } else if (side!=CblasLeft && side!=CblasRight) {
        info = 2;
    } else if (upLo!=CblasLower && upLo!=CblasUpper) {
        info = 3;
    } else if (transA!=CblasNoTrans && transA!=CblasTrans
            && transA!=AtlasConj && transA!=CblasConjTrans)
    {
        info = 4;
    } else if (diag!=CblasUnit && diag!=CblasNonUnit) {
        info = 5;
    } else {
        if (order==CblasColMajor) {
            if (m<0) {
                info = 6;
            } else if (n<0) {
                info = 7;
            } else if (ldA<std::max(1,numRowsA)) {
                info = 10;
            } else if (ldB<std::max(1,m)) {
                info = 12;
            }
        } else {
            if (n<0) {
                info = 6;
            } else if (m<0) {
                info = 7;
            } else if (ldA<std::max(1,numRowsA)) {
                info = 10;
            } else if (ldB<std::max(1,n)) {
                info = 12;
            }
        }
    }

    if (info!=0) {
        CBLAS(xerbla)(info, "cblas_dtrmm", "");
    }

    if (order==CblasColMajor) {
        ULMBLAS(dtrmm)(side, upLo, transA, diag, m, n, alpha, A, ldA, B, ldB);
    } else {
        side   = (side==CblasLeft)  ? CblasRight : CblasLeft;
        upLo   = (upLo==CblasUpper) ? CblasLower : CblasUpper;
        ULMBLAS(dtrmm)(side, upLo, transA, diag, n, m, alpha, A, ldA, B, ldB);
    }
}

} // extern "C"
