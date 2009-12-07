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

#ifndef FORTRAN_NAME
# define FORTRAN_NAME(lower,upper) lower##_
#endif
typedef long fortran_int;
typedef long int fortran_logical;
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
MatrixOpStatus Copy(const SubVectorView<TrivialReadableVectorView<TVector<double_complex,TAlloc> > > &src, TVector<double_complex,TAlloc> &dst){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zcopy,ZCOPY)(src.size(), src.Raw(), 1, dst.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Copy(const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &src, TVector<double_complex,TAlloc> &dst){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zcopy,ZCOPY)(src.size(), src.Raw(), 1, dst.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Copy(const TVector<double_complex,TAlloc> &src, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &dst){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zcopy,ZCOPY)(src.size(), src.Raw(), 1, dst.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Copy(const SubVectorView<TrivialReadableVectorView<TVector<double_complex,TAlloc> > > &src, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &dst){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zcopy,ZCOPY)(src.size(), src.Raw(), 1, dst.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Copy(const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &src, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &dst){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zcopy,ZCOPY)(src.size(), src.Raw(), 1, dst.Raw(), 1);
	return OK;
}

template <class TAlloc>
MatrixOpStatus Copy(const SubVectorView<ColumnView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > > &src, const ColumnView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zcopy,ZCOPY)(dst.size(), src.Raw(), 1, dst.Raw(), 1);
	return OK;
}

template <class TAlloc>
MatrixOpStatus Copy(const TVector<double_complex,TAlloc> &src, DiagonalView<TrivialWritableMatrixView<TDiagonalMatrix<double_complex,TAlloc> > > &dst){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zcopy,ZCOPY)(src.size(), src.Raw(), 1, &dst.GetMutable(0), 1);
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
MatrixOpStatus Copy(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &src, TMatrix<double_complex,TAlloc> &dst){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < dst.Cols(); ++j){
		FORTRAN_NAME(zcopy,ZCOPY)(dst.Rows(), src.Raw()+j*src.LeadingDimension(), 1, &(dst.GetMutable(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Copy(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &src, TMatrix<double_complex,TAlloc> &dst){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < dst.Cols(); ++j){
		FORTRAN_NAME(zcopy,ZCOPY)(dst.Rows(), src.Raw()+j*src.LeadingDimension(), 1, &(dst.GetMutable(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Copy(const TMatrix<double_complex,TAlloc> &src, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < dst.Cols(); ++j){
		FORTRAN_NAME(zcopy,ZCOPY)(dst.Rows(), src.Raw()+j*src.LeadingDimension(), 1, &(dst.GetMutable(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Copy(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &src, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < dst.Cols(); ++j){
		FORTRAN_NAME(zcopy,ZCOPY)(dst.Rows(), src.Raw()+j*src.LeadingDimension(), 1, &(dst.GetMutable(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Copy(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &src, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < dst.Cols(); ++j){
		FORTRAN_NAME(zcopy,ZCOPY)(dst.Rows(), src.Raw()+j*src.LeadingDimension(), 1, &(dst.GetMutable(0,j)), 1);
	}
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
		std::uninitialized_fill_n(&(dst.GetMutable(0,j)), dst.Rows(), value);
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
		dst.Set(i, value);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Fill(const DiagonalView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &D, const double_complex &value){
	for(size_t i = 0; i < D.size(); ++i){
		D.Set(i, value);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Fill(const DiagonalView<SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > > &D, const double_complex &value){
	for(size_t i = 0; i < D.size(); ++i){
		D.Set(i, value);
	}
	return OK;
}




extern "C" void FORTRAN_NAME(zswap,ZSWAP)(const fortran_int &N, double_complex *x, const fortran_int &incx,
               double_complex *y, const fortran_int &incy );
template <class TAlloc>
MatrixOpStatus TransposeInPlace(TMatrix<double_complex,TAlloc> &A){
	assert(A.Rows() == A.Cols());
	const size_t n = A.Cols()-1;
	for(size_t j = 0; j < n; ++j){
		FORTRAN_NAME(zswap,ZSWAP)(n-j, &(A(j+1,j)), 1, &(A(j,j+1)), A.LeadingDimension());
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
MatrixOpStatus Scale(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const double &scale){
	for(size_t j = 0; j < A.Cols(); ++j){
		FORTRAN_NAME(zdscal,ZDSCAL)(A.Rows(), scale, &A.GetMutable(0,j), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Scale(TMatrix<double_complex,TAlloc> &A, const double_complex &scale){
	FORTRAN_NAME(zscal,ZSCAL)(A.Rows()*A.Cols(), scale, A.Raw(), 1);
	return OK;
}
template <class TAlloc>
MatrixOpStatus Scale(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const double_complex &scale){
	for(size_t j = 0; j < A.Cols(); ++j){
		FORTRAN_NAME(zscal,ZSCAL)(A.Rows(), scale, &A.GetMutable(0,j), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Scale(TVector<double_complex,TAlloc> &v, const double_complex &scale){
	for(size_t j = 0; j < v.size(); ++j){
		FORTRAN_NAME(zscal,ZSCAL)(v.size(), scale, v.Raw(), 1);
	}
	return OK;
}












extern "C" double FORTRAN_NAME(dznrm2,DZNRM2)(const fortran_int &N, const double_complex *x, const fortran_int &incx);
template <class TAlloc>
double Norm2(const TVector<double_complex,TAlloc> &v){
	return FORTRAN_NAME(dznrm2,DZNRM2)(v.size(), v.Raw(), 1);
}
template <class TAlloc>
double Norm2(const ColumnView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &v){
	return FORTRAN_NAME(dznrm2,DZNRM2)(v.size(), &v.Get(0), 1);
}
template <class TAlloc>
double FrobeniusNorm(const TMatrix<double_complex,TAlloc> &A){
	return FORTRAN_NAME(dznrm2,DZNRM2)(A.Rows()*A.Cols(), A.Raw(), 1);
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
		FORTRAN_NAME(zscal,ZSCAL)(A.Rows(), D[j], &(A.GetMutable(0,j)), 1);
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
		FORTRAN_NAME(zscal,ZSCAL)(A.Cols(), D[i], &(A.GetMutable(i,0)), A.LeadingDimension());
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TDiagonalMatrix<double_complex> &D, TVector<double_complex,TAlloc> &X){
	assert(D.size() == X.size());
	for(size_t i = 0; i < X.size(); ++i){
		X[i] *= D[i];
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Mult(const TDiagonalMatrix<double_complex> &D, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &X){
	assert(D.size() == X.size());
	for(size_t i = 0; i < X.size(); ++i){
		X.Set(i, D[i] * X.Get(i));
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
MatrixOpStatus Add(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &src, TMatrix<double_complex,TAlloc> &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < src.Cols(); ++j){
		FORTRAN_NAME(zaxpy,ZAXPY)(src.Rows(), scale_src, src.Raw()+j*src.LeadingDimension(), 1, &(dst.GetMutable(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &src, TMatrix<double_complex,TAlloc> &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < src.Cols(); ++j){
		FORTRAN_NAME(zaxpy,ZAXPY)(src.Rows(), scale_src, src.Raw()+j*src.LeadingDimension(), 1, &(dst.GetMutable(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const TMatrix<double_complex,TAlloc> &src, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < src.Cols(); ++j){
		FORTRAN_NAME(zaxpy,ZAXPY)(src.Rows(), scale_src, &(src(0,j)), 1, &(dst.GetMutable(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &src, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < src.Cols(); ++j){
		FORTRAN_NAME(zaxpy,ZAXPY)(src.Rows(), scale_src, src.Raw()+j*src.LeadingDimension(), 1, &(dst.GetMutable(0,j)), 1);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &src, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	for(size_t j = 0; j < src.Cols(); ++j){
		FORTRAN_NAME(zaxpy,ZAXPY)(src.Rows(), scale_src, src.Raw()+j*src.LeadingDimension(), 1, &(dst.GetMutable(0,j)), 1);
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
MatrixOpStatus Add(const SubVectorView<TrivialReadableVectorView<TVector<double_complex,TAlloc> > > &src, TVector<double_complex,TAlloc> &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zaxpy,ZAXPY)(src.size(), scale_src, src.Raw(), src.Stride(), dst.Raw(), 1);
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
MatrixOpStatus Add(const SubVectorView<TrivialReadableVectorView<TVector<double_complex,TAlloc> > > &src, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zaxpy,ZAXPY)(src.size(), scale_src, src.Raw(), src.Stride(), dst.Raw(), dst.Stride());
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &src, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	FORTRAN_NAME(zaxpy,ZAXPY)(src.size(), scale_src, src.Raw(), src.Stride(), dst.Raw(), dst.Stride());
	return OK;
}

#ifdef USING_TONES_VECTOR
template <class TAlloc>
MatrixOpStatus Add(const TOnesVector<double_complex> &src, const DiagonalView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	for(size_t i = 0; i < dst.size(); ++i){
		dst.Set(i, dst.Get(i) + scale_src);
	}
	return OK;
}
#endif
template <class TAlloc>
MatrixOpStatus Add(const TDiagonalMatrix<double_complex> &src, const DiagonalView<SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	for(size_t i = 0; i < dst.size(); ++i){
		dst.Set(i, dst.Get(i) + scale_src*src[i]);
	}
	return OK;
}
template <class TAlloc>
MatrixOpStatus Add(const DiagonalView<TrivialWritableMatrixView<TDiagonalMatrix<double_complex> > > &src, const DiagonalView<SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > > &dst, const double_complex &scale_src = double_complex(1)){
	assert(src.size() == dst.size());
	for(size_t i = 0; i < dst.size(); ++i){
		dst.Set(i, dst.Get(i) + scale_src*src.Get(i));
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
MatrixOpStatus Mult(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const TVector<double_complex,TAlloc> &X, TVector<double_complex,TAlloc> &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
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
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const SubVectorView<TrivialReadableVectorView<TVector<double_complex,TAlloc> > > &X, TVector<double_complex,TAlloc> &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), X.Stride(),
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
MatrixOpStatus Mult(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubVectorView<TrivialReadableVectorView<TVector<double_complex,TAlloc> > > &X, TVector<double_complex,TAlloc> &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
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
MatrixOpStatus Mult(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const TVector<double_complex,TAlloc> &X, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
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
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const SubVectorView<TrivialReadableVectorView<TVector<double_complex,TAlloc> > > &X, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
	FORTRAN_NAME(zgemv,ZGEMV)("N", A.Rows(), A.Cols(),
		scale_AX, A.Raw(), A.LeadingDimension(),
		X.Raw(), X.Stride(),
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
MatrixOpStatus Mult(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubVectorView<TrivialReadableVectorView<TVector<double_complex,TAlloc> > > &X, const SubVectorView<TrivialWritableVectorView<TVector<double_complex,TAlloc> > > &Y, const double_complex &scale_AX = double_complex(1), const double_complex &scale_Y = double_complex(0)){
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
MatrixOpStatus Mult(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const TMatrix<double_complex,TAlloc> &B, TMatrix<double_complex,TAlloc> &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
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
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &B, TMatrix<double_complex,TAlloc> &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
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
MatrixOpStatus Mult(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &B, TMatrix<double_complex,TAlloc> &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
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
MatrixOpStatus Mult(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const TMatrix<double_complex,TAlloc> &B, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
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
MatrixOpStatus Mult(const TMatrix<double_complex,TAlloc> &A, const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &B, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
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
MatrixOpStatus Mult(const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &A, const SubMatrixView<TrivialReadableMatrixView<TMatrix<double_complex,TAlloc> > > &B, const SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &C, const double_complex &scale_AB = double_complex(1), const double_complex &scale_C = double_complex(0)){
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
template <class TAlloc>
MatrixOpStatus InvertDestructive(SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, TMatrix<double_complex,TAlloc> &Ainv){
	Fill(Ainv, double_complex(0));
	Fill(Diagonal(Ainv), double_complex(1));
	return SolveDestructive(A, Ainv);
}
template <class TAlloc>
MatrixOpStatus InvertDestructive(TMatrix<double_complex,TAlloc> &A, SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &Ainv){
	Fill(Ainv, double_complex(0));
	Fill(Diagonal(Ainv), double_complex(1));
	return SolveDestructive(A, Ainv);
}
template <class TAlloc>
MatrixOpStatus InvertDestructive(SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &A, SubMatrixView<TrivialWritableMatrixView<TMatrix<double_complex,TAlloc> > > &Ainv){
	Fill(Ainv, double_complex(0));
	Fill(Diagonal(Ainv), double_complex(1));
	return SolveDestructive(A, Ainv);
}



// inv(A + UCV) = inv(A) - inv(A)*U*inv(inv(C) + V*inv(A)*U)*V*inv(A)
template <class TAlloc>
MatrixOpStatus UpdateInverse(TMatrix<double_complex,TAlloc> &Ainv, const TMatrix<double_complex,TAlloc> &U, const TMatrix<double_complex,TAlloc> &C, const TMatrix<double_complex,TAlloc> &V){
	const size_t k = C.Rows();
	const size_t n = Ainv.Rows();
	assert(Ainv.Rows() == U.Rows());
	assert(U.Cols() == C.Rows());
	assert(C.Cols() == V.Rows());
	assert(V.Cols() == Ainv.Cols());
	assert(C.Rows() == C.Cols());
	assert(Ainv.Rows() == Ainv.Cols());
	
	TMatrix<double_complex,TAlloc> iAU(n,k);
	TMatrix<double_complex,TAlloc> ViA(k,n), iC(k,k);
	Copy(C, SubMatrix(ViA, 0,0,k,k));
	InvertDestructive(ViA, iC); // inv(C) in iC
	
	Mult(Ainv, U, iAU);
	Mult(V, iAU, iC, double_complex(1), double_complex(1)); // iC = inv(C) + V*inv(A)*U
	MatrixOpStatus ret = InvertDestructive(iC, SubMatrix(ViA, 0,0,k,k));
	Copy(SubMatrix(ViA, 0,0,k,k), iC); // ic = inv(inv(C) + V*inv(A)*U)
	
	Mult(V, Ainv, ViA);
	TMatrix<double_complex,TAlloc> iAUB(n,k);
	Mult(iAU, iC, iAUB);
	Mult(iAUB, ViA, Ainv, complex_t(-1), complex_t(1));
	return ret;
}








extern "C" void FORTRAN_NAME(zgeev,ZGEEV)( const char *jobvl, const char *jobvr, const fortran_int &N,
           double_complex *a, const fortran_int &lda, double_complex *w,
           double_complex *vl, const fortran_int &ldvl, double_complex *vr,
           const fortran_int &ldvr, double_complex *work, const fortran_int &lwork,
           double *rwork, fortran_int &info );
template <class TAlloc>
inline MatrixOpStatus Eigensystem(TMatrix<double_complex,TAlloc> &A, TVector<double_complex,TAlloc> &Eval, TMatrix<double_complex,TAlloc> &EvecLeft, TMatrix<double_complex,TAlloc> &EvecRight){
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == EvecRight.Rows());
	assert(A.Cols() == EvecRight.Cols());
	assert(A.Rows() == Eval.size());
	assert(A.Rows() == EvecLeft.Rows());
	assert(A.Cols() == EvecLeft.Cols());

	TAlloc allocator;
	fortran_int info(1);
	fortran_int lwork(1+2*(int)A.Rows());
	typedef typename TAlloc::template rebind<double>::other double_allocator;
	double_allocator dallocator(allocator);
	double *rwork = dallocator.allocate(2*A.Rows());
	double_complex *work = allocator.allocate(lwork);
	FORTRAN_NAME(zgeev,ZGEEV)("V", "V", A.Rows(), A.Raw(), A.LeadingDimension(), Eval.Raw(), EvecLeft.Raw(), EvecLeft.LeadingDimension(), EvecRight.Raw(), EvecRight.LeadingDimension(), work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, 2*A.Rows());
	return (0 == info) ? OK : UNKNOWN_ERROR;
}
template <class TAlloc>
inline MatrixOpStatus Eigensystem(TMatrix<double_complex,TAlloc> &A, TVector<double_complex,TAlloc> &Eval, TMatrix<double_complex,TAlloc> &Evec){
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == Evec.Rows());
	assert(A.Cols() == Evec.Cols());
	assert(A.Rows() == Eval.size());

	TAlloc allocator;
	fortran_int info(1);
	fortran_int lwork(1+2*(int)A.Rows());
	typedef typename TAlloc::template rebind<double>::other double_allocator;
	double_allocator dallocator(allocator);
	double *rwork = dallocator.allocate(2*A.Rows());
	double_complex *work = allocator.allocate(lwork);
	FORTRAN_NAME(zgeev,ZGEEV)("N", "V", A.Rows(), A.Raw(), A.LeadingDimension(), Eval.Raw(), NULL, 1, Evec.Raw(), Evec.LeadingDimension(), work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, 2*A.Rows());
	return (0 == info) ? OK : UNKNOWN_ERROR;
}
template <class TAlloc>
inline MatrixOpStatus Eigensystem(TMatrix<double_complex,TAlloc> &A, TVector<double_complex,TAlloc> &Eval){
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == Eval.size());

	TAlloc allocator;
	fortran_int info(1);
	fortran_int lwork(1+2*(int)A.Rows());
	typedef typename TAlloc::template rebind<double>::other double_allocator;
	double_allocator dallocator(allocator);
	double *rwork = dallocator.allocate(2*A.Rows());
	double_complex *work = allocator.allocate(lwork);
	FORTRAN_NAME(zgeev,ZGEEV)("N", "N", A.Rows(), A.Raw(), A.LeadingDimension(), Eval.Raw(), NULL, 1, NULL, 1, work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, 2*A.Rows());
	return (0 == info) ? OK : UNKNOWN_ERROR;
}


extern "C" void FORTRAN_NAME(zggev,ZGGEV)( const char *jobvl, const char *jobvr, const fortran_int &N,
           double_complex *a, const fortran_int &lda, double_complex *b, const fortran_int &ldb,
           double_complex *alpha, double_complex *beta,
           double_complex *vl, const fortran_int &ldvl, double_complex *vr,
           const fortran_int &ldvr, double_complex *work, const fortran_int &lwork,
           double *rwork, fortran_int &info );
extern "C" void FORTRAN_NAME(zggevx,ZGGEVX)(
	const char *balanc, const char *jobvl, const char *jobvr, const char *sense,
	const fortran_int &n, double_complex *a, const fortran_int &lda,
	double_complex *b, const fortran_int &ldb,
	double_complex *alpha, double_complex *beta, 
	double_complex *vl, const fortran_int &ldvl,
	double_complex *vr, const fortran_int &ldvr, 
	fortran_int &ilo, fortran_int &ihi,
	double *lscale, double *rscale, 
	double &abnrm, double &bbnrm, double *rconde, double *rcondv,
	double_complex *work, const fortran_int &lwork, double *rwork, 
	fortran_int *iwork, fortran_logical *bwork, fortran_int &info);
template <class TAlloc>
inline MatrixOpStatus GeneralizedEigensystem(TMatrix<double_complex,TAlloc> &A, TMatrix<double_complex,TAlloc> &B, TVector<double_complex,TAlloc> &alpha, TVector<double_complex,TAlloc> &beta, TMatrix<double_complex,TAlloc> &EvecLeft, TMatrix<double_complex,TAlloc> &EvecRight){
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == EvecRight.Rows());
	assert(A.Cols() == EvecRight.Cols());
	assert(A.Rows() == alpha.size());
	assert(A.Rows() == beta.size());
	assert(A.Rows() == EvecLeft.Rows());
	assert(A.Cols() == EvecLeft.Cols());

	TAlloc allocator;
	fortran_int info(1);
	fortran_int lwork(2*(int)A.Rows());
	typedef typename TAlloc::template rebind<double>::other double_allocator;
	double_allocator dallocator(allocator);
	double *rwork = dallocator.allocate(8*A.Rows());
	double_complex *work = allocator.allocate(lwork);
	FORTRAN_NAME(zggev,ZGGEV)("V", "V", A.Rows(), A.Raw(), A.LeadingDimension(), B.Raw(), B.LeadingDimension(), alpha.Raw(), beta.Raw(), EvecLeft.Raw(), EvecLeft.LeadingDimension(), EvecRight.Raw(), EvecRight.LeadingDimension(), work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, 8*A.Rows());
	return (0 == info) ? OK : UNKNOWN_ERROR;
}
template <class TAlloc>
inline MatrixOpStatus GeneralizedEigensystem(TMatrix<double_complex,TAlloc> &A, TMatrix<double_complex,TAlloc> &B, TVector<double_complex,TAlloc> &alpha, TVector<double_complex,TAlloc> &beta, TMatrix<double_complex,TAlloc> &Evec){
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == Evec.Rows());
	assert(A.Cols() == Evec.Cols());
	assert(A.Rows() == alpha.size());
	assert(A.Rows() == beta.size());

	TAlloc allocator;
	typedef typename TAlloc::template rebind<double>::other DAlloc;
	DAlloc dallocator(allocator);

	fortran_int info(1);
	
	fortran_int lwork(2*(int)A.Rows());
	double *rwork = dallocator.allocate(8*A.Rows());
	double_complex *work = allocator.allocate(lwork);
	FORTRAN_NAME(zggev,ZGGEV)("N", "V", A.Rows(), A.Raw(), A.LeadingDimension(), B.Raw(), B.LeadingDimension(), alpha.Raw(), beta.Raw(), NULL, 1, Evec.Raw(), Evec.LeadingDimension(), work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, 8*A.Rows());
	return (0 == info) ? OK : UNKNOWN_ERROR;
}
template <class TAlloc>
inline MatrixOpStatus GeneralizedEigensystem(TMatrix<double_complex,TAlloc> &A, TMatrix<double_complex,TAlloc> &B, TVector<double_complex,TAlloc> &alpha, TVector<double_complex,TAlloc> &beta, TMatrix<double_complex,TAlloc> &Evec, double *EvalTol){
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == Evec.Rows());
	assert(A.Cols() == Evec.Cols());
	assert(A.Rows() == alpha.size());
	assert(A.Rows() == beta.size());

	TAlloc allocator;
	typedef typename TAlloc::template rebind<double>::other DAlloc;
	typedef typename TAlloc::template rebind<fortran_int>::other IAlloc;
	typedef typename TAlloc::template rebind<fortran_logical>::other LAlloc;
	DAlloc dallocator(allocator);
	IAlloc iallocator(allocator);
	LAlloc lallocator(allocator);

	fortran_int info(1);
	/*
	fortran_int lwork(2*(int)A.Rows());
	double *rwork = dallocator.allocate(8*A.Rows());
	double_complex *work = allocator.allocate(lwork);
	FORTRAN_NAME(zggev,ZGGEV)("N", "V", A.Rows(), A.Raw(), A.LeadingDimension(), B.Raw(), B.LeadingDimension(), alpha.Raw(), beta.Raw(), NULL, 1, Evec.Raw(), Evec.LeadingDimension(), work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, 8*A.Rows());
	*/
	fortran_int lwork(2*(int)A.Rows() * ((NULL == EvalTol) ? 1:2));
	fortran_int lrwork(6*(int)A.Rows());
	fortran_int liwork((int)A.Rows() + 2);
	fortran_int lbwork((int)A.Rows());
	double_complex *work = allocator.allocate(lwork);
	double *rwork = dallocator.allocate(lrwork);
	fortran_int *iwork = NULL;
	fortran_logical *bwork = NULL;
	if(NULL != EvalTol){
		bwork = lallocator.allocate(lbwork);
	}else{
		iwork = iallocator.allocate(liwork);
	}

	fortran_int ilo, ihi;
	double abnrm, bbnrm;
	double *lscale = dallocator.allocate((int)A.Rows());
	double *rscale = dallocator.allocate((int)A.Rows());
	FORTRAN_NAME(zggevx,ZGGEVX)("S", "N", "V", ((NULL == EvalTol) ? "N" : "E"),
		A.Rows(), A.Raw(), A.LeadingDimension(), B.Raw(), B.LeadingDimension(), alpha.Raw(), beta.Raw(),
		NULL, 1, Evec.Raw(), Evec.LeadingDimension(),
		ilo, ihi, lscale, rscale, abnrm, bbnrm, EvalTol, NULL,
		work, lwork, rwork, iwork, bwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, lrwork);
	dallocator.deallocate(lscale, A.Rows());
	dallocator.deallocate(rscale, A.Rows());
	if(NULL != EvalTol){
		lallocator.deallocate(bwork, lbwork);
		double scale = std::sqrt(abnrm*abnrm + bbnrm*bbnrm) * std::numeric_limits<double>::epsilon();
		for(size_t i = 0; i < A.Rows(); ++i){
			EvalTol[i] = scale/EvalTol[i];
		}
	}else{
		iallocator.deallocate(iwork, liwork);
	}
	return (0 == info) ? OK : UNKNOWN_ERROR;
}
template <class TAlloc>
inline MatrixOpStatus GeneralizedEigensystem(TMatrix<double_complex,TAlloc> &A, TMatrix<double_complex,TAlloc> &B, TVector<double_complex,TAlloc> &alpha, TVector<double_complex,TAlloc> &beta){
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == alpha.size());
	assert(A.Rows() == beta.size());

	TAlloc allocator;
	fortran_int info(1);
	fortran_int lwork(2*(int)A.Rows());
	typedef typename TAlloc::template rebind<double>::other double_allocator;
	double_allocator dallocator(allocator);
	double *rwork = dallocator.allocate(8*A.Rows());
	double_complex *work = allocator.allocate(lwork);
	FORTRAN_NAME(zggev,ZGGEV)("N", "N", A.Rows(), A.Raw(), A.LeadingDimension(), B.Raw(), B.LeadingDimension(), alpha.Raw(), beta.Raw(), NULL, 1, NULL, 1, work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, 8*A.Rows());
	return (0 == info) ? OK : UNKNOWN_ERROR;
}

extern "C" void FORTRAN_NAME(zgesvd,ZGESVD)(const char *jobU, const char *jobV, const fortran_int &M, const fortran_int &N,
                   const double_complex *A, const fortran_int &lda,
                   double *S,
                   const double_complex *U, const fortran_int &ldU,
                   const double_complex *VT, const fortran_int &ldVT,
                   double_complex *work, const fortran_int &lwork, double *rwork, fortran_int &info);
extern "C" void FORTRAN_NAME(zgesdd,ZGESDD)(const char *jobU, const char *jobV, const fortran_int &M, const fortran_int &N,
                   const double_complex *A, const fortran_int &lda,
                   double *S,
                   const double_complex *U, const fortran_int &ldU,
                   const double_complex *VT, const fortran_int &ldVT,
                   double_complex *work, const fortran_int &lwork, double *rwork, fortran_int *iwork, fortran_int &info);
template <class TAlloc>
inline MatrixOpStatus SingularValueDecompositionDestructive(const TMatrix<double_complex,TAlloc> &A, TVector<double,typename TAlloc::template rebind<double>::other> &sigma){
	const size_t min_dim = (A.Rows() < A.Cols()) ? A.Rows() : A.Cols();
	const size_t max_dim = (A.Rows() > A.Cols()) ? A.Rows() : A.Cols();
	fortran_int info(1);
	fortran_int lwork(2*min_dim+max_dim);
	fortran_int lrwork = 5*min_dim;
	TAlloc allocator;
	double_complex *work = allocator.allocate(lwork);
	
	typedef typename TAlloc::template rebind<double>::other double_allocator;
	double_allocator dallocator(allocator);
	double *rwork = dallocator.allocate(lrwork);
	sigma.Resize(min_dim);
	
	FORTRAN_NAME(zgesvd,ZGESVD)("N", "N", A.Rows(), A.Cols(), A.Raw(), A.LeadingDimension(), sigma.Raw(), NULL, 1, NULL, 1, work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, lrwork);
	return (info == 0) ? OK : UNKNOWN_ERROR;
}
template <class TAlloc>
inline MatrixOpStatus SingularValueDecompositionDestructive(const TMatrix<double_complex,TAlloc> &A, TMatrix<double_complex,TAlloc> &U, TVector<double,typename TAlloc::template rebind<double>::other> &sigma, TMatrix<double_complex,TAlloc> &Vt){
	const size_t min_dim = (A.Rows() < A.Cols()) ? A.Rows() : A.Cols();
	const size_t max_dim = (A.Rows() > A.Cols()) ? A.Rows() : A.Cols();
	fortran_int info(1);
	fortran_int lwork(2*min_dim+max_dim);
	fortran_int lrwork = 5*min_dim;
	TAlloc allocator;
	double_complex *work = allocator.allocate(lwork);
	
	typedef typename TAlloc::template rebind<double>::other double_allocator;
	double_allocator dallocator(allocator);
	double *rwork = dallocator.allocate(lrwork);
	sigma.Resize(min_dim);
	U.Resize(A.Rows(), min_dim);
	Vt.Resize(A.Cols(), min_dim);
	
	FORTRAN_NAME(zgesvd,ZGESVD)("A", "A", A.Rows(), A.Cols(), A.Raw(), A.LeadingDimension(), sigma.Raw(), U.Raw(), U.LeadingDimension(), Vt.Raw(), Vt.LeadingDimension(), work, lwork, rwork, info);
	allocator.deallocate(work, lwork);
	dallocator.deallocate(rwork, lrwork);
	return (info == 0) ? OK : UNKNOWN_ERROR;
}
template <class TAlloc>
inline double Norm2(TMatrix<double_complex,TAlloc> &A){
	TVector<double,typename TAlloc::template rebind<double>::other> sigma;
	TMatrix<double_complex,TAlloc> Acopy(A);
	SingularValueDecompositionDestructive(Acopy, sigma);
	return sigma[0];
}
template <class TAlloc>
inline double ConditionNumber(TMatrix<double_complex,TAlloc> &A){
	TVector<double,typename TAlloc::template rebind<double>::other> sigma;
	TMatrix<double_complex,TAlloc> Acopy(A);
	SingularValueDecompositionDestructive(Acopy, sigma);
	return sigma[0]/sigma[sigma.size()-1];
}

//// UnitaryProcrustes(A) - Replaces A with the nearest (Frobenius norm) unitary matrix, A'*A = I
//   Solved by taking SVD of A, and replacing the singular values with 1's
template <class TAlloc>
inline MatrixOpStatus UnitaryProcrustes(TMatrix<double_complex,TAlloc> &A){
	TMatrix<double_complex,TAlloc> U, Vt;
	TVector<double_complex,TAlloc> sigma;
	MatrixOpStatus ret = SingularValueDecompositionDestructive(A, U, sigma, Vt);
	Mult(U, Vt, A);
	
	return ret;
}


extern "C" void FORTRAN_NAME(zpotrf,ZPOTRF)(const char *uplo, const fortran_int &N,
                   double_complex *A, const fortran_int &lda,
                   fortran_int &info);
extern "C" void FORTRAN_NAME(ztrmm,ZTRMM)(const char *side, const char *uplo, const char *transA, const char *diag,
                   const fortran_int &M, const fortran_int &N, const double_complex &alpha,
                   double_complex *A, const fortran_int &ldA,
                   double_complex *B, const fortran_int &ldB);
extern "C" void FORTRAN_NAME(ztrtri,ZTRTRI)(const char *uplo, const char *diag,
                   const fortran_int &N, double_complex *A, const fortran_int &ldA,
                   fortran_int &info);
//// GeneralizedProcrustes(A, B, C) - Replaces A with nearest (Frobenius norm) matrix such that A'*B*A = C
//   Solved with Cholesky factorization of B and C (both must be Hermitian positive definite)
//   Let B = L*L' and C = M*M', and iL and iM are inverses of L and M, respectively
//   Then we need
//     A'*U'*U*A = V'*V
//     iV'*A'*U'*U*A*iV = I
//   Let Q = U*A*iV and solve UnitaryProcrustes(Q)
//   Then A is iU*Q*V
template <class TAlloc>
inline MatrixOpStatus GeneralizedUnitaryProcrustes(TMatrix<std::complex<double>,TAlloc> &A, const TMatrix<std::complex<double>,TAlloc> &B, const TMatrix<std::complex<double>,TAlloc> &C){
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
	
	FORTRAN_NAME(zpotrf,ZPOTRF)("U", U.Rows(), U.LeadingDimension(), info);
	FORTRAN_NAME(zpotrf,ZPOTRF)("U", V.Rows(), V.LeadingDimension(), info);
	
	FORTRAN_NAME(ztrmm,ZTRMM)("L", "U", "N", "N", A.Rows(), A.Cols(), complex_t(1), U.Raw(), U.LeadingDimension(), A.Raw(), A.LeadingDimension()); // A <- U*A;
	
	TMatrix<std::complex<double> > iV(V);
	FORTRAN_NAME(ztrtri,ZTRTRI)("U", "N", iV.Rows(), iV.LeadingDimension(), info);
	
	FORTRAN_NAME(ztrmm,ZTRMM)("R", "U", "N", "N", A.Rows(), A.Cols(), complex_t(1), iV.Raw(), iV.LeadingDimension(), A.Raw(), A.LeadingDimension()); // A is now U*A*iV where the last A is the original A
	
	MatrixOpStatus ret = UnitaryProcrustes(A);
	
	FORTRAN_NAME(ztrtri,ZTRTRI)("U", "N", U.Rows(), U.LeadingDimension(), info); // U is now iU
	
	FORTRAN_NAME(ztrmm,ZTRMM)("L", "U", "N", "N", A.Rows(), A.Cols(), complex_t(1), U.Raw(), U.LeadingDimension(), A.Raw(), A.LeadingDimension());
	FORTRAN_NAME(ztrmm,ZTRMM)("R", "U", "N", "N", A.Rows(), A.Cols(), complex_t(1), V.Raw(), V.LeadingDimension(), A.Raw(), A.LeadingDimension());
	
	return ret;
}

}; // namespace MatrixOpsFortran

#endif // _MATRIX_OPS_FORTRAN_H_

