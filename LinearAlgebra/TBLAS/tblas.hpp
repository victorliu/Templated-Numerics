#ifndef _TBLAS_H_
#define _TBLAS_H_

#include <complex>

// In the naive implementations, level N functions are only allowed to call
// level N-1 or lower functions.

namespace TBLAS{

// A general template interface
// Specializations should typedef real_t
template <typename T>
class TypeTraits{
public:
	// typedef ... real_t;
	// virtual static T Conjugate(T x)
};

template <>
class TypeTraits<double>{
public:
	typedef double real_t;
	static double Conjugate(double x){ return x; }
};

template <>
class TypeTraits<std::complex<double> >{
public:
	typedef double real_t;
	static std::complex<double> Conjugate(std::complex<double> x){ return std::conj(x); }
};

/*
class Op{
private:
	static const char cNone = 'N';
	static const char cTranspose = 'T';
	static const char cConjugateTranspose = 'C';
public:
	static const char *None(){ return &cNone; }
	static const char *Transpose(){ return &cTranspose; }
	static const char *ConjugateTranspose(){ return &cConjugateTranspose; }
};
*/

template <class T>
class OpBase{
public:
	static const T None;
	static const T Transpose;
	static const T ConjugateTranspose;
}
template <>
class OpBase<char>{
	static const char None = 'N';
	static const char Transpose = 'T';
	static const char ConjugateTranspose = 'C';
};
typedef OpBase<char> Op;

class Side{
private:
	static const char cLeft = 'L';
	static const char cRight = 'R';
public:
	static const char *Left(){ return &cLeft; }
	static const char *Right(){ return &cRight; }
};

class UpLo{
private:
	static const char cUpper = 'U';
	static const char cLower = 'L';
public:
	static const char *Upper(){ return &cUpper; }
	static const char *Lower(){ return &cLower; }
};

class Diag{
private:
	static const char cUnit = 'U';
	static const char cNonUnit = 'N';
public:
	static const char *Unit(){ return &cUnit; }
	static const char *NonUnit(){ return &cNonUnit; }
};

//////////////////////////////////////////////
//////////////// LEVEL 1 BLAS ////////////////
//////////////////////////////////////////////

template <typename ScalarT>
void TBLAS_NAME(rotg,ROTG)(
  ScalarT &a,
  ScalarT &b,
  ScalarT &c,
  ScalarT &s);
#include <blas/tblas_1_rotg.hpp>

template <typename ScalarRealT>
void TBLAS_NAME(rotmg,ROTMG)(
  ScalarRealT &d1,
  ScalarRealT &d2,
  ScalarRealT &x1,
  const ScalarRealT &y1,
  ScalarRealT param[5]);
#include <blas/tblas_1_rotmg.hpp>

template <typename ScalarRealT>
void TBLAS_NAME(rot,ROT)(
  const size_t &n,
  ScalarRealT *x, const size_t &incx,
  ScalarRealT *y, const size_t &incy,
  const ScalarRealT &c, const ScalarRealT &s);
#include <blas/tblas_1_rot.hpp>

template <typename ScalarRealT>
void TBLAS_NAME(rotm,ROTM)(
  const size_t &n,
  ScalarRealT *x, const size_t &incx,
  ScalarRealT *y, const size_t &incy,
  const ScalarRealT param[5]);
#include <blas/tblas_1_rotm.hpp>

template <typename ScalarT>
void TBLAS_NAME(swap,SWAP)(
  const size_t &n,
  ScalarT *x, const size_t &incx,
  ScalarT *y, const size_t &incy);
#include <blas/tblas_1_swap.hpp>

template <typename ScalarT>
void TBLAS_NAME(scal,SCAL)(
  const size_t &n,
  const ScalarT &alpha,
  ScalarT *x, const size_t &incx);
#include <blas/tblas_1_scal.hpp>

template <typename ScalarT>
void TBLAS_NAME(copy,COPY)(
  const size_t &n,
  const ScalarT *src, const size_t &incsrc,
  ScalarT *dst, const size_t &incdst);
#include <blas/tblas_1_copy.hpp>

template <typename ScalarT>
void TBLAS_NAME(axpy,AXPY)(
  const size_t &n,
  const ScalarT &alpha,
  const ScalarT *x, const size_t &incx,
  ScalarT *y, const size_t &incy);
#include <blas/tblas_1_axpy.hpp>

#ifdef TBLAS_EXTENSIONS

template <typename ScalarT>
void TBLAS_NAME(axpby,AXPBY)( // y <- alpha*x + beta*y
  const size_t &n,
  const ScalarT &alpha,
  const ScalarT *x, const size_t &incx,
  const ScalarT &beta,
  ScalarT *y, const size_t &incy);
#include <blas/tblas_1x_axpby.hpp>

template <typename DiagT, typename ScalarT>
void TBLAS_NAME(discal,DISCAL)( // x <- diag(D)*x
  const size_t &n,
  const DiagT *D, const size_t &incd,
  ScalarT *x, const size_t &incx);
#include <blas/tblas_1x_discal.hpp>

template <typename DiagT, typename ScalarT>
void TBLAS_NAME(dxpby,DXPBY)( // y <- diag(D)*x + beta*y
  const size_t &n,
  const DiagT *D, const size_t &incd,
  const ScalarT *x, const size_t &incx,
  const ScalarT &beta,
  ScalarT *y, const size_t &incy);
#include <blas/tblas_1x_dxpby.hpp>

#endif // TBLAS_EXTENSIONS

template <typename ScalarT1, typename ScalarT2>
ScalarT1 TBLAS_NAME(dot,DOT)( // note: return type is type of left argument
  const size_t &n,
  const ScalarT1 *x, const size_t &incx,
  const ScalarT2 *y, const size_t &incy);
#include <blas/tblas_1_dot.hpp>

template <typename ScalarT1, typename ScalarT2>
ScalarT1 TBLAS_NAME(dotu,DOTU)( // same as dot
  const size_t &n,
  const ScalarT1 *x, const size_t &incx,
  const ScalarT2 *y, const size_t &incy);
#include <blas/tblas_1_dotu.hpp>

template <typename ScalarT1, typename ScalarT2>
ScalarT1 TBLAS_NAME(dotc,DOTC)( // note: return type is type of left argument
  const size_t &n,
  const ScalarT1 *x, const size_t &incx,
  const ScalarT2 *y, const size_t &incy);
#include <blas/tblas_1_dotc.hpp>

template <typename ScalarT>
ScalarT TBLAS_NAME(sdsdot,SDSDOT)(
  const size_t &n,
  const ScalarT *x, const size_t &incx,
  const ScalarT *y, const size_t &incy);
#include <blas/tblas_1_sdsdot.hpp>

template <typename ScalarT>
typename TBLAS_TRAITS(ScalarT)::real_t TBLAS_NAME(nrm2,NRM2)(
  const size_t &n,
  ScalarT *x, const size_t &incx);
#include <blas/tblas_1_nrm2.hpp>

template <typename ScalarT>
typename TBLAS_TRAITS(ScalarT)::real_t TBLAS_NAME(asum,ASUM)(
  const size_t &n,
  const ScalarT *x, const size_t &incx);
#include <blas/tblas_1_asum.hpp>

template <typename ScalarT>
size_t TBLAS_NAME(iamax,IAMAX)(
  const size_t &n,
  const ScalarT *x, const size_t &incx);
#include <blas/tblas_1_iamax.hpp>

//////////////////////////////////////////////
//////////////// LEVEL 2 BLAS ////////////////
//////////////////////////////////////////////

// size of A is ALWAYS m by n (even when trans is 'T')
template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(gemv,GEMV)(
  const char *trans,
  const size_t &m, const size_t &n,
  const TypeAlpha& alpha,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx,
  const TypeBeta& beta,
  const TypeY *y, const size_t &incy);
#include <blas/tblas_2_gemv.hpp>

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(gbmv,GBMV)(
  const char *trans,
  const size_t &m, const size_t &n,
  const size_t &kl, const size_t &ku,
  const TypeAlpha& alpha,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx,
  const TypeBeta& beta,
  const TypeY *y, const size_t &incy);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(hemv,HEMV)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha& alpha,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx,
  const TypeBeta& beta,
  const TypeY *y, const size_t &incy);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(hbmv,HBMV)(
  const char *uplo,
  const size_t &n,
  const size_t &k,
  const TypeAlpha& alpha,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx,
  const TypeBeta& beta,
  const TypeY *y, const size_t &incy);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(hpmv,HPMV)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha& alpha,
  const TypeA *AP,
  const TypeX *x, const size_t &incx,
  const TypeBeta& beta,
  const TypeY *y, const size_t &incy);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(symv,SYMV)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha& alpha,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx,
  const TypeBeta& beta,
  const TypeY *y, const size_t &incy);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(sbmv,SBMV)(
  const char *uplo,
  const size_t &n,
  const size_t &k,
  const TypeAlpha& alpha,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx,
  const TypeBeta& beta,
  const TypeY *y, const size_t &incy);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(spmv,SPMV)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha& alpha,
  const TypeA *AP,
  const TypeX *x, const size_t &incx,
  const TypeBeta& beta,
  const TypeY *y, const size_t &incy);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(trmv,TRMV)(
  const char *uplo, const char *trans, const char *diag,
  const size_t &n,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(tbmv,TBMV)(
  const char *uplo, const char *trans, const char *diag,
  const size_t &n,
  const size_t &k,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(tpmv,TPMV)(
  const char *uplo, const char *trans, const char *diag,
  const size_t &n,
  const TypeA *AP,
  const TypeX *x, const size_t &incx);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(trsv,TRSV)(
  const char *uplo, const char *trans, const char *diag,
  const size_t &n,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(tbsv,TBSV)(
  const char *uplo, const char *trans, const char *diag,
  const size_t &n,
  const size_t &k,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx);

template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(tpsv,TPSV)(
  const char *uplo, const char *trans, const char *diag,
  const size_t &n,
  const TypeA *AP,
  const TypeX *x, const size_t &incx);


#ifdef TBLAS_EXTENSIONS

// A <- diag(D)*A or A <- A*diag(D)
template <typename TypeD, typename TypeA>
void TBLAS_NAME(gedm,GEDM)(
  const char *side,
  const size_t &m, const size_t &n,
  const TypeD *D, const size_t &incd,
  const TypeA *A, const size_t &lda);

#endif // TBLAS_EXTENSIONS


template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(ger,GER)(
  const size_t &m, const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeX *y, const size_t &incy,
  const TypeA *A, const size_t &lda);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(geru,GERU)(
  const size_t &m, const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeX *y, const size_t &incy,
  const TypeA *A, const size_t &lda);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(gerc,GERC)(
  const size_t &m, const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeX *y, const size_t &incy,
  const TypeA *A, const size_t &lda);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(her,HER)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeA *A, const size_t &lda);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(hpr,HPR)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeA *AP);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(her2,HER2)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeX *y, const size_t &incy,
  const TypeA *A, const size_t &lda);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(hpr2,HPR2)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeX *y, const size_t &incy,
  const TypeA *AP);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(syr,SYR)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeA *A, const size_t &lda);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(spr,SPR)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeA *AP);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(syr2,SYR2)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeX *y, const size_t &incy,
  const TypeA *A, const size_t &lda);

template <typename TypeAlpha, typename TypeX, typename TypeY, typename TypeA>
void TBLAS_NAME(spr2,SPR2)(
  const char *uplo,
  const size_t &n,
  const TypeAlpha &alpha,
  const TypeX *x, const size_t &incx,
  const TypeX *y, const size_t &incy,
  const TypeA *AP);



//////////////////////////////////////////////
//////////////// LEVEL 3 BLAS ////////////////
//////////////////////////////////////////////


template <typename TypeAlpha, typename TypeA, typename TypeB,
          typename TypeBeta, typename TypeC>
void TBLAS_NAME(gemm,GEMM)(
  const char *transa, const char *transb,
  const size_t &m, const size_t &n, const size_t &k,
  const TypeAlpha &alpha,
  const TypeA *A, const size_t &lda,
  const TypeB *B, const size_t &ldb,
  const TypeBeta &beta,
  const TypeC *C, const size_t &ldc);

template <typename TypeAlpha, typename TypeA, typename TypeB,
          typename TypeBeta, typename TypeC>
void TBLAS_NAME(symm,SYMM)(
  const char *side, const char *uplo,
  const size_t &m, const size_t &n,
  const TypeAlpha &alpha,
  const TypeA *A, const size_t &lda,
  const TypeB *B, const size_t &ldb,
  const TypeBeta &beta,
  const TypeC *C, const size_t &ldc);

template <typename TypeAlpha, typename TypeA, typename TypeB,
          typename TypeBeta, typename TypeC>
void TBLAS_NAME(hemm,HEMM)(
  const char *side, const char *uplo,
  const size_t &m, const size_t &n,
  const TypeAlpha &alpha,
  const TypeA *A, const size_t &lda,
  const TypeB *B, const size_t &ldb,
  const TypeBeta &beta,
  const TypeC *C, const size_t &ldc);

template <typename TypeAlpha, typename TypeA,
          typename TypeBeta, typename TypeC>
void TBLAS_NAME(syrk,SYRK)(
  const char *uplo, const char *trans,
  const size_t &n, const size_t &k,
  const TypeAlpha &alpha,
  const TypeA *A, const size_t &lda,
  const TypeBeta &beta,
  const TypeC *C, const size_t &ldc);

template <typename TypeAlpha, typename TypeA,
          typename TypeBeta, typename TypeC>
void TBLAS_NAME(herk,HERK)(
  const char *uplo, const char *trans,
  const size_t &n, const size_t &k,
  const TypeAlpha &alpha,
  const TypeA *A, const size_t &lda,
  const TypeBeta &beta,
  const TypeC *C, const size_t &ldc);

template <typename TypeAlpha, typename TypeA, typename TypeB,
          typename TypeBeta, typename TypeC>
void TBLAS_NAME(syr2k,SYR2K)(
  const char *uplo, const char *trans,
  const size_t &n, const size_t &k,
  const TypeAlpha &alpha,
  const TypeA *A, const size_t &lda,
  const TypeB *B, const size_t &ldb,
  const TypeBeta &beta,
  const TypeC *C, const size_t &ldc);

template <typename TypeAlpha, typename TypeA, typename TypeB,
          typename TypeBeta, typename TypeC>
void TBLAS_NAME(her2k,HER2K)(
  const char *uplo, const char *trans,
  const size_t &n, const size_t &k,
  const TypeAlpha &alpha,
  const TypeA *A, const size_t &lda,
  const TypeB *B, const size_t &ldb,
  const TypeBeta &beta,
  const TypeC *C, const size_t &ldc);

template <typename TypeAlpha, typename TypeA, typename TypeB>
void TBLAS_NAME(trmm,TRMM)(
  const char *side, const char *uplo, const char *trans, const char *diag,
  const size_t &m, const size_t &n,
  const TypeAlpha &alpha,
  const TypeA *A, const size_t &lda,
  const TypeB *B, const size_t &ldb);

}; // namespace TBLAS

#endif // _TBLAS_H_
