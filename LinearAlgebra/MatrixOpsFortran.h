#ifndef _MATRIX_OPS_FORTRAN_H_
#define _MATRIX_OPS_FORTRAN_H_

#include <memory>

#ifndef MATRIX_OP_STATUS_DEFINED
#define MATRIX_OP_STATUS_DEFINED
enum MatrixOpStatus{
	OK,
	DIMENSION_MISMATCH,
	SINGULAR_MATRIX,
	UNKNOWN_ERROR
};
#endif

namespace MatrixOpsFortran{

#define FORTRAN_NAME(lower,upper) lower##_
typedef long fortran_int;
typedef std::complex<double> double_complex;

extern "C" void FORTRAN_NAME(zcopy,ZCOPY)(
	const fortran_int &N,
	const double_complex *x, const fortran_int &incx,
	double_complex *y, const fortran_int &incy
);
template <class TAlloc>
MatrixOpStatus Copy(const TVector<double_complex,TAlloc> &src, TVector<double_complex,TAlloc> &dst){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zcopy,ZCOPY)(src.size(), src.Raw(), 1, dst.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Copy(const TMatrix<double_complex,TAlloc> &src, TMatrix<double_complex,TAlloc> &dst){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	FORTRAN_NAME(zcopy,ZCOPY)(src.Rows()*src.Cols(), src.Raw(), 1, dst.Raw(), 1);
	return OK;
}



template <class TAlloc>
MatrixOpStatus Fill(TMatrix<double_complex,TAlloc> &dst, const double_complex &value){
	std::uninitialized_fill_n(dst.Raw(), dst.Rows()*dst.Cols(), value);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Fill(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst, const double_complex &value){
	for(size_t j = 0; j < dst.Cols(); ++j){
		std::uninitialized_fill_n(&(dst(0,j)), dst.Rows(), value);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Fill(TVector<double_complex,TAlloc> &dst, const double_complex &value){
	std::uninitialized_fill_n(dst.Raw(), dst.size(), value);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Fill(const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &dst, const double_complex &value){
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] = value;
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Fill(const DiagonalView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &D, const double_complex &value){
	for(size_t i = 0; i < D.size(); ++i){
		D[i] = value;
	}
	return OK;
}




extern "C" void FORTRAN_NAME(zdscal,ZDSCAL)( const fortran_int &N, const double &alpha,
                double_complex *x, const fortran_int &incx );
extern "C" void FORTRAN_NAME(zscal,ZSCAL)( const fortran_int &N, const double_complex &alpha,
               double_complex *x, const fortran_int &incx );
template <class TAlloc>
MatrixOpStatus Scale(TMatrix<double_complex,TAlloc> &A, const double &scale){
	FORTRAN_NAME(zdscal,ZDSCAL)(A.Rows()*A.Cols(), scale, A.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Scale(TMatrix<double_complex,TAlloc> &A, const double_complex &scale){
	FORTRAN_NAME(zscal,ZSCAL)(A.Rows()*A.Cols(), scale, A.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(TMatrix<double_complex,TAlloc> &A, const TDiagonalMatrix<double_complex> &D){
	assert(D.size() == A.Cols());
	for(size_t j = 0; j < A.Cols(); ++j){
		FORTRAN_NAME(zscal,ZSCAL)(A.Rows(), D[j], &(A(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const TDiagonalMatrix<double_complex> &D){
	assert(D.size() == A.Cols());
	for(size_t j = 0; j < A.Cols(); ++j){
		FORTRAN_NAME(zscal,ZSCAL)(A.Rows(), D[j], &(A(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TDiagonalMatrix<double_complex> &D, TMatrix<double_complex,TAlloc> &A){
	assert(D.size() == A.Cols());
	for(size_t i = 0; i < A.Rows(); ++i){
		FORTRAN_NAME(zscal,ZSCAL)(A.Cols(), D[i], &(A(i,0)), A.LeadingDimension());
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TDiagonalMatrix<double_complex> &D, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A){
	assert(D.size() == A.Cols());
	for(size_t i = 0; i < A.Rows(); ++i){
		FORTRAN_NAME(zscal,ZSCAL)(A.Cols(), D[i], &(A(i,0)), A.LeadingDimension());
	}
	return OK;
}








extern "C" void FORTRAN_NAME(zaxpy,ZAXPY)(
	const fortran_int &N, const double_complex &alpha,
	const double_complex *x, const fortran_int &incx,
	double_complex *y, const fortran_int &incy );
template <class TAlloc>
MatrixOpStatus Add(const TMatrix<double_complex,TAlloc> &src, TMatrix<double_complex,TAlloc> &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	FORTRAN_NAME(zaxpy,ZAXPY)(src.Rows()*src.Cols(), scale_src, src.Raw(), 1, dst.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &src, TMatrix<double_complex,TAlloc> &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < src.Cols(); ++j){
		FORTRAN_NAME(zaxpy,ZAXPY)(src.Rows(), scale_src, &(src(0,j)), 1, &(dst(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const TMatrix<double_complex,TAlloc> &src, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < src.Cols(); ++j){
		FORTRAN_NAME(zaxpy,ZAXPY)(src.Rows(), scale_src, &(src(0,j)), 1, &(dst(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &src, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < src.Cols(); ++j){
		FORTRAN_NAME(zaxpy,ZAXPY)(src.Rows(), scale_src, &(src(0,j)), 1, &(dst(0,j)), 1);
	}
	return OK;
}

template <class TAlloc>
MatrixOpStatus Add(const TVector<double_complex,TAlloc> &src, TVector<double_complex,TAlloc> &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zaxpy,ZAXPY)(src.size(), scale_src, src.Raw(), 1, dst.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &src, TVector<double_complex,TAlloc> &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zaxpy,ZAXPY)(src.size(), scale_src, src.Raw(), src.Stride(), dst.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const TVector<double_complex,TAlloc> &src, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zaxpy,ZAXPY)(src.size(), scale_src, src.Raw(), 1, dst.Raw(), dst.Stride());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &src, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zaxpy,ZAXPY)(src.size(), scale_src, src.Raw(), src.Stride(), dst.Raw(), dst.Stride());
	return OK;
}

template <class TAlloc>
MatrixOpStatus Add(const TOnesVector<double_complex> &src, const DiagonalView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] += scale_src;
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const TDiagonalMatrix<double_complex> &src, const DiagonalView<SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] += scale_src*src[i];
	}
	return OK;
}






extern "C" void FORTRAN_NAME(zgemv,ZGEMV)( const char *trans, const fortran_int &M, const fortran_int &N,
           const double_complex &alpha, const double_complex *A, const fortran_int &lda,
           const double_complex *X, const fortran_int &incX,
           const double_complex &beta, double_complex *Y, const fortran_int &incY);
template <class TAlloc>
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const TVector<double_complex,TAlloc> &X, TVector<double_complex,TAlloc> &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), 1,
		scale_Y, Y.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const TVector<double_complex,TAlloc> &X, TVector<double_complex,TAlloc> &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), 1,
		scale_Y, Y.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &X, TVector<double_complex,TAlloc> &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), X.Stride(),
		scale_Y, Y.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &X, TVector<double_complex,TAlloc> &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), X.Stride(),
		scale_Y, Y.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const TVector<double_complex,TAlloc> &X, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), 1,
		scale_Y, Y.Raw(), Y.Stride());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const TVector<double_complex,TAlloc> &X, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), 1,
		scale_Y, Y.Raw(), Y.Stride());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &X, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), X.Stride(),
		scale_Y, Y.Raw(), Y.Stride());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &X, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), X.Stride(),
		scale_Y, Y.Raw(), Y.Stride());
	return OK;
}






extern "C" void FORTRAN_NAME(zgemm,ZGEMM)( const char *transA, const char *transB, const fortran_int &M, const fortran_int &N, const fortran_int &K,
           const double_complex &alpha, const double_complex *A, const fortran_int &lda,
           const double_complex *B, const fortran_int &ldb,
           const double_complex &beta, double_complex *C, const fortran_int &ldc);
template <class TAlloc>
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const TMatrix<double_complex,TAlloc> &B, TMatrix<double_complex,TAlloc> &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
	FORTRAN_NAME(zgemm,ZGEMM)("N", "N", A.Rows(), B.Cols(), A.Cols(),
		scale_AB, A.Raw(), A.LeadingDimension(),
		B.Raw(), B.LeadingDimension(),
		scale_C, C.Raw(), C.LeadingDimension());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const TMatrix<double_complex,TAlloc> &B, TMatrix<double_complex,TAlloc> &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
	FORTRAN_NAME(zgemm,ZGEMM)("N", "N", A.Rows(), B.Cols(), A.Cols(),
		scale_AB, A.Raw(), A.LeadingDimension(),
		B.Raw(), B.LeadingDimension(),
		scale_C, C.Raw(), C.LeadingDimension());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &B, TMatrix<double_complex,TAlloc> &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
	FORTRAN_NAME(zgemm,ZGEMM)("N", "N", A.Rows(), B.Cols(), A.Cols(),
		scale_AB, A.Raw(), A.LeadingDimension(),
		B.Raw(), B.LeadingDimension(),
		scale_C, C.Raw(), C.LeadingDimension());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &B, TMatrix<double_complex,TAlloc> &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
	FORTRAN_NAME(zgemm,ZGEMM)("N", "N", A.Rows(), B.Cols(), A.Cols(),
		scale_AB, A.Raw(), A.LeadingDimension(),
		B.Raw(), B.LeadingDimension(),
		scale_C, C.Raw(), C.LeadingDimension());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const TMatrix<double_complex,TAlloc> &B, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
	FORTRAN_NAME(zgemm,ZGEMM)("N", "N", A.Rows(), B.Cols(), A.Cols(),
		scale_AB, A.Raw(), A.LeadingDimension(),
		B.Raw(), B.LeadingDimension(),
		scale_C, C.Raw(), C.LeadingDimension());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const TMatrix<double_complex,TAlloc> &B, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
	FORTRAN_NAME(zgemm,ZGEMM)("N", "N", A.Rows(), B.Cols(), A.Cols(),
		scale_AB, A.Raw(), A.LeadingDimension(),
		B.Raw(), B.LeadingDimension(),
		scale_C, C.Raw(), C.LeadingDimension());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &B, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
	FORTRAN_NAME(zgemm,ZGEMM)("N", "N", A.Rows(), B.Cols(), A.Cols(),
		scale_AB, A.Raw(), A.LeadingDimension(),
		B.Raw(), B.LeadingDimension(),
		scale_C, C.Raw(), C.LeadingDimension());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &B, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
	FORTRAN_NAME(zgemm,ZGEMM)("N", "N", A.Rows(), B.Cols(), A.Cols(),
		scale_AB, A.Raw(), A.LeadingDimension(),
		B.Raw(), B.LeadingDimension(),
		scale_C, C.Raw(), C.LeadingDimension());
	return OK;
}








extern "C" void FORTRAN_NAME(zgesv,ZGESV)( const fortran_int &M, const fortran_int &NRHS,
           double_complex *A, const fortran_int &lda,
           fortran_int *ipiv,
           double_complex *B, const fortran_int &ldb,
           fortran_int &info);
template <class TAlloc>
MatrixOpStatus SolveDestructive(TMatrix<double_complex,TAlloc> &A, TMatrix<double_complex,TAlloc> &X){
	fortran_int *ipiv = new fortran_int[A.Rows()];
	fortran_int info = 1;
	FORTRAN_NAME(zgesv,ZGESV)(A.Rows(), X.Cols(), A.Raw(), A.LeadingDimension(), ipiv, X.Raw(), X.LeadingDimension(), info);
	delete [] ipiv;
	if(0 == info){ return OK; }
	else{ return SINGULAR_MATRIX; }
}
template <class TAlloc>
MatrixOpStatus SolveDestructive(SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, TMatrix<double_complex,TAlloc> &X){
	fortran_int *ipiv = new fortran_int[A.Rows()];
	fortran_int info = 1;
	FORTRAN_NAME(zgesv,ZGESV)(A.Rows(), X.Cols(), A.Raw(), A.LeadingDimension(), ipiv, X.Raw(), X.LeadingDimension(), info);
	delete [] ipiv;
	if(0 == info){ return OK; }
	else{ return SINGULAR_MATRIX; }
}
template <class TAlloc>
MatrixOpStatus SolveDestructive(TMatrix<double_complex,TAlloc> &A, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &X){
	fortran_int *ipiv = new fortran_int[A.Rows()];
	fortran_int info = 1;
	FORTRAN_NAME(zgesv,ZGESV)(A.Rows(), X.Cols(), A.Raw(), A.LeadingDimension(), ipiv, X.Raw(), X.LeadingDimension(), info);
	delete [] ipiv;
	if(0 == info){ return OK; }
	else{ return SINGULAR_MATRIX; }
}
template <class TAlloc>
MatrixOpStatus SolveDestructive(SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &X){
	fortran_int *ipiv = new fortran_int[A.Rows()];
	fortran_int info = 1;
	FORTRAN_NAME(zgesv,ZGESV)(A.Rows(), X.Cols(), A.Raw(), A.LeadingDimension(), ipiv, X.Raw(), X.LeadingDimension(), info);
	delete [] ipiv;
	if(0 == info){ return OK; }
	else{ return SINGULAR_MATRIX; }
}



template <class TAlloc>
MatrixOpStatus SolveDestructive(TMatrix<double_complex,TAlloc> &A, TVector<double_complex,TAlloc> &X){
	fortran_int *ipiv = new fortran_int[A.Rows()];
	fortran_int info = 1;
	FORTRAN_NAME(zgesv,ZGESV)(A.Rows(), 1, A.Raw(), A.LeadingDimension(), ipiv, X.Raw(), X.size(), info);
	delete [] ipiv;
	if(0 == info){ return OK; }
	else{ return SINGULAR_MATRIX; }
}
template <class TAlloc>
MatrixOpStatus SolveDestructive(SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, TVector<double_complex,TAlloc> &X){
	fortran_int *ipiv = new fortran_int[A.Rows()];
	fortran_int info = 1;
	FORTRAN_NAME(zgesv,ZGESV)(A.Rows(), 1, A.Raw(), A.LeadingDimension(), ipiv, X.Raw(), X.size(), info);
	delete [] ipiv;
	if(0 == info){ return OK; }
	else{ return SINGULAR_MATRIX; }
}
template <class TAlloc>
MatrixOpStatus SolveDestructive(TMatrix<double_complex,TAlloc> &A, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &X){
	fortran_int *ipiv = new fortran_int[A.Rows()];
	fortran_int info = 1;
	FORTRAN_NAME(zgesv,ZGESV)(A.Rows(), 1, A.Raw(), A.LeadingDimension(), ipiv, X.Raw(), X.size(), info);
	delete [] ipiv;
	if(0 == info){ return OK; }
	else{ return SINGULAR_MATRIX; }
}
template <class TAlloc>
MatrixOpStatus SolveDestructive(SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &X){
	fortran_int *ipiv = new fortran_int[A.Rows()];
	fortran_int info = 1;
	FORTRAN_NAME(zgesv,ZGESV)(A.Rows(), 1, A.Raw(), A.LeadingDimension(), ipiv, X.Raw(), X.size(), info);
	delete [] ipiv;
	if(0 == info){ return OK; }
	else{ return SINGULAR_MATRIX; }
}







template <class TAlloc>
MatrixOpStatus InvertDestructive(TMatrix<double_complex,TAlloc> &A, TMatrix<double_complex,TAlloc> &Ainv){
	Fill(Ainv, double_complex(0));
	Fill(Diagonal(Ainv), double_complex(1));
	return SolveDestructive(A, Ainv);
}










extern "C" void FORTRAN_NAME(zgeev,ZGEEV)( const char *jobvl, const char *jobvr, const fortran_int &N,
           double_complex *a, const fortran_int &lda, double_complex *w,
           double_complex *vl, const fortran_int &ldvl, double_complex *vr,
           const fortran_int &ldvr, double_complex *work, const fortran_int &lwork,
           double *rwork, fortran_int &info );
template <class TAlloc>
inline MatrixOpStatus Eigensystem(TMatrix<double_complex,TAlloc> &A, TVector<double_complex,TAlloc> &Eval, TMatrix<double_complex,TAlloc> &Evec){
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == Evec.Rows());
	assert(A.Cols() == Evec.Cols());
	assert(A.Rows() == Eval.size());

	TAlloc allocator;
	fortran_int info(1);
	fortran_int lwork(1+2*(int)A.Rows());
	double *rwork = TAlloc::rebind<double>::other(allocator).allocate(2*A.Rows());
	double_complex *work = allocator.allocate(lwork);
	FORTRAN_NAME(zgeev,ZGEEV)("N", "V", A.Rows(), A.Raw(), A.LeadingDimension(), Eval.Raw(), NULL, 1, Evec.Raw(), Evec.LeadingDimension(), work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	TAlloc::rebind<double>::other(allocator).deallocate(rwork, 2*A.Rows());
	return (0 == info) ? OK : UNKNOWN_ERROR;
}


/*
#ifdef USE_COMPLEX_MATRICESS
extern "C" zgesvd_(const char *jobU, const char *jobV, const long &M, const long &N,
                   const std::complex<double> *A, const long &lda,
                   std::complex<double> *S,
                   const std::complex<double> *U, const long &ldU,
                   const std::complex<double> *VT, const long &ldVT,
                   std::complex<double> *work, const long &lwork, long &info);
//// UnitaryProcrustes(A) - Replaces A with the nearest (Frobenius norm) unitary matrix, A'*A = I
//   Solved by taking SVD of A, and replacing the singular values with 1's
template <class TAlloc>
inline MatrixOpStatus UnitaryProcrustes(TMatrix<std::complex<double>,TAlloc> &A){
	assert(A.Rows() == A.Cols());
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != A.Cols()){ return DIMENSION_MISMATCH; }
#endif
	long info(1);
	long lwork(3*(int)A.Rows());
	double *rwork = new double[5*A.Rows()];
	TMatrix<std::complex<double> > U(A.Rows(), A.Cols());
	TMatrix<std::complex<double> > VT(A.Cols(), A.Rows());
	TVector<std::complex<double> > S(A.Rows());
	std::complex<double> *work = new std::complex<double>[lwork];
	zgesvd_("A", "A", A.Rows(), A.Cols(), A.Raw(), A.LeadingDimension(),
		S.Raw(), U.Raw(), U.LeadingDimension(),
		VT.Raw(), VT.LeadingDimension(),
		work, lwork, rwork, info);
	delete [] work;
	delete [] rwork;
	
	Mult(U,VT,A);
	
	return (0 == info) ? OK : UNKNOWN_ERROR;
}


extern "C" zpotrf_(const char *uplo, const long &N,
                   std::complex<double> *A, const long &lda,
                   long &info);
extern "C" ztrmm_(const char *side, const char *uplo, const char *transA, const char *diag,
                   const long &M, const long &N, const std::complex<double> &alpha,
                   std::complex<double> *A, const long &ldA,
                   std::complex<double> *B, const long &ldB);
extern "C" ztrtri_(const char *uplo, const char *diag,
                   const long &N, std::complex<double> *A, const long &ldA,
                   long &info);
//// GeneralizedProcrustes(A, B, C) - Replaces A with nearest (Frobenius norm) matrix such that A'*B*A = C
//   Solved with Cholesky factorization of B and C (both must be Hermitian positive definite)
//   Let B = L*L' and C = M*M', and iL and iM are inverses of L and M, respectively
//   Then we need
//     A'*U'*U*A = V'*V
//     iV'*A'*U'*U*A*iV = I
//   Let Q = U*A*iV and solve UnitaryProcrustes(Q)
//   Then A is iU*Q*V
template <class TAlloc>
inline MatrixOpStatus GeneralizedProcrustes(TMatrix<std::complex<double>,TAlloc> &A, const TMatrix<std::complex<double>,TAlloc> &B, const TMatrix<std::complex<double>,TAlloc> &C){
	assert(A.Rows() == A.Cols());
	assert(B.Rows() == B.Cols());
	assert(C.Rows() == C.Cols());
	assert(A.Rows() == B.Rows());
	assert(A.Rows() == C.Rows());
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != A.Cols() || B.Rows() != B.Cols() || C.Rows() != C.Cols() || A.Rows() != B.Rows() || A.Rows() != C.Rows()){ return DIMENSION_MISMATCH; }
#endif
	TMatrix<std::complex<double> > U(B), V(C), temp(A.Rows(), A.Cols());
	long info;
	
	zpotrf_("U", U.Rows(), U.LeadingDimension(), info);
	zpotrf_("U", V.Rows(), V.LeadingDimension(), info);
	
	ztrmm_("L", "U", "N", "N", A.Rows(), A.Cols(), complex_t(1), U.Raw(), U.LeadingDimension(), A.Raw(), A.LeadingDimension()); // A <- U*A;
	
	TMatrix<std::complex<double> > iV(V);
	ztrtri_("U", "N", iV.Rows(), iV.LeadingDimension(), info);
	
	ztrmm_("R", "U", "N", "N", A.Rows(), A.Cols(), complex_t(1), iV.Raw(), iV.LeadingDimension(), A.Raw(), A.LeadingDimension()); // A is now U*A*iV where the last A is the original A
	
	MatrixOpStatus ret = UnitaryProcrustes(A);
	
	ztrtri_("U", "N", U.Rows(), U.LeadingDimension(), info); // U is now iU
	
	ztrmm_("L", "U", "N", "N", A.Rows(), A.Cols(), complex_t(1), iU.Raw(), iU.LeadingDimension(), A.Raw(), A.LeadingDimension());
	ztrmm_("R", "U", "N", "N", A.Rows(), A.Cols(), complex_t(1), V.Raw(), V.LeadingDimension(), A.Raw(), A.LeadingDimension());
	
	return ret;
}
#endif // USE_COMPLEX_MATRICESS
*/

}; // namespace MatrixOpsFortran

#endif // _MATRIX_OPS_FORTRAN_H_

