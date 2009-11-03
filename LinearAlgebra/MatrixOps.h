#ifndef _MATRIX_OPS_H_
#define _MATRIX_OPS_H_

// Preprocessor flags:
//   USE_MATRIX_OPS_ASSERTS - Turn on assertion checking
//   USE_MATRIX_OPS_TBLAS   - Use optimized BLAS routines when possible

#ifdef USE_MATRIX_OPS_ASSERTS
# include <cassert>
#endif

namespace MatrixOps{

enum MatrixOpStatus{
	OK,
	DIMENSION_MISMATCH,
	SINGULAR_MATRIX,
	UNKNOWN_ERROR
};

// Fundamental operations
//   Copy(src, dst)
//   Fill(dst, value)
//   Scale(X, scaleX)
//   Add(X, YPlusX, scaleX = 1)
//   Mult(A, B, CPlusATimesB, scaleATimesB = 1, scale_C = 0)

//// Fill
template <class T>
MatrixOpStatus Fill(TMatrixBase<T> &dst, const T &value){
	for(size_t i = 0; i < dst.Rows(); ++i){
		for(size_t j = 0; j < dst.Cols(); ++j){
			dst(i,j) = value;
		}
	}
	return OK;
}
template <class T>
MatrixOpStatus Fill(TVectorBase<T> &dst, const T &value){
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] = value;
	}
	return OK;
}

//// Mult

template <class TA, class TB, class TC>
MatrixOpStatus Mult(
	const TA &A, const TB &B, TC &C,
	typename TC::value_type &scale_AB = typename TC::value_type(1),
	typename TC::value_type &scale_C = typename TC::value_type(0)
){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
#endif
	for(size_t i = 0; i < A.Rows(); ++i){
		for(size_t j = 0; j < B.Cols(); ++j){
			TC::value_type sum(0);
			for(size_t k = 0; k < A.Cols(); ++k){
				sum += A(i,k)*B(k,j);
			}
			C(i,j) = scale_AB * sum + scale_C * C(i,j);
		}
	}
	return OK;
}

#ifdef USE_MATRIX_OPS_TBLAS

template <>
MatrixOpStatus Mult<
	TMatrix<double>,
	TMatrix<double>,
	TMatrix<double>
>(
	const TMatrix<double> &A, const TMatrix<double> &B, TMatrix<double> &C,
	double &scale_AB = 1.0,
	double &scale_C = 0.0
){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
#endif
	TBLAS::TBLAS_NAME(gemm,GEMM)(
		TBLAS::Op::None, TBLAS::Op::None,
		A.Rows(), B.Cols(), A.Cols(),
		
	return OK;
}

template <>
MatrixOpStatus Mult<
	TransposeView<TMatrixView<double> >,
	TMatrix<double>,
	TMatrix<double>
>(
	const TransposeView<TMatrixView<double> > &A, const TMatrix<double> &B, TMatrix<double> &C,
	double &scale_AB = 1.0,
	double &scale_C = 0.0
){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
#endif
	TBLAS::TBLAS_NAME(gemm,GEMM)(
	return OK;
}

#endif // USE_MATRIX_OPS_TBLAS

}; // namespace MatrixOps

#endif // _MATRIX_OPS_H_
