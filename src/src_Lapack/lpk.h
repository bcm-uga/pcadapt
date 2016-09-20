#ifndef LPK_H
#define LPK_H

#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"

int dgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c, integer *ldc);

int dgemmlig_(char *transa, char *transb, integer *m, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c, integer *ldc);

int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *info);

int dgesdd_(char *jobz, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *iwork, integer *info);

int dgetrf_(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info);

int dgetri_(integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);

int dorgtr_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

int dstedc_(char *compz, integer *n, doublereal *d__, doublereal *e, doublereal *z__, integer *ldz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

int dsvdc_(doublereal *x, integer *ldx, integer *n, integer *p, doublereal *s, doublereal *e, doublereal *u, integer *ldu, doublereal *v, integer *ldv, doublereal *work, integer *job, integer *info);

int dsytrd_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *d__, doublereal *e, doublereal *tau, doublereal *work, integer *lwork, integer *info);

int dsyevr_(char *jobz, char *range, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

#endif //LPK_H
