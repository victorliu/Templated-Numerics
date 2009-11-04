#ifndef _MATRIX_OPS_H_
#define _MATRIX_OPS_H_

// Preprocessor flags:
//   USE_MATRIX_OPS_CHECKS  - If not turned on, always returns OK
//   USE_MATRIX_ASSERTS - Turn on assertion checking
//   USE_MATRIX_OPS_TBLAS   - Use optimized BLAS routines when possible
//   USE_COMPLEX_MATRICES

// These are fully generic matrix operations. Particular routines for
// different classes are contained in their respective headers.

#ifdef USE_MATRIX_ASSERTS
# include <cassert>
#endif

#ifdef USE_COMPLEX_MATRICES
# include <complex>
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
//   Dot(X, Y)
//   ConjugateDot(Xc, Y)
//   Scale(X, scaleX)
//   Add(X, YPlusX, scaleX = 1)
//   Rank1Update(A, X, Yt, scaleXY = 1)
//   Rank1Update(A, X, scaleXX = 1)
//   Rank2Update(A, X, Yt, scaleXY = 1)
//   Mult(A, B, CPlusATimesB, scaleATimesB = 1, scale_C = 0)
//
// Sophisticated operations
//   Solve(A, B, X) - A*X == B
//   SolveLeastSquares(A, B, X) - A*X == B
//   Eigenvalues(A, L)
//   Eigensystem(A, L, X) - A*X == X*diag(L)
//   SingularValueDecomposition(A, U, S, V) - A == U*S*V
//   LUDecomposition(A, L, U) - A == L*U
//   CholeskyDecomposition(A, L) - A == L*L'

//// Copy

template <class TSRC, class TDST>
MatrixOpStatus Copy(const TMatrixBase<TSRC> &src, TMatrixBase<TDST> &dst){
#ifdef USE_MATRIX_OPS_CHECKS
	if(src.Rows() != dst.Rows() && src.Cols() != dst.Cols()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(src.Rows() != dst.Rows());
	assert(src.Cols() != dst.Cols());
#endif
	for(size_t j = 0; j < dst.Cols(); ++j){
		for(size_t i = 0; i < dst.Rows(); ++i){
			dst(i,j) = TDST(src(i,j));
		}
	}
	return OK;
}
template <class TSRC, class TDST>
MatrixOpStatus Copy(const TVectorBase<TSRC> &src, TVectorBase<TDST> &dst){
#ifdef USE_MATRIX_OPS_CHECKS
	if(src.size() != dst.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(src.size() != dst.size());
#endif
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] = TDST(src[i]);
	}
	return OK;
}

//// Fill

template <class T, class V>
MatrixOpStatus Fill(TMatrixBase<T> &dst, const V &value){
	for(size_t j = 0; j < dst.Cols(); ++j){
		for(size_t i = 0; i < dst.Rows(); ++i){
			dst(i,j) = typename TMatrixBase<T>::value_type(value);
		}
	}
	return OK;
}
template <class T, class V>
MatrixOpStatus Fill(TVectorBase<T> &dst, const V &value){
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] = typename TVectorBase<T>::value_type(value);
	}
	return OK;
}

//// Dot

template <class TX, class TY>
TX Dot(const TVectorBase<TX> &x, const TVectorBase<TY> &y){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(x.size() != y.size());
#endif
	TX sum(0);
	for(size_t i = 0; i < x.size(); ++i){
		sum += x[i]*y[i];
	}
	return sum;
}

//// ConjugateDot

#ifdef USE_COMPLEX_MATRICES

template <class TXC, class TY>
TXC ConjugateDot(const TVectorBase<TXC> &xc, const TVectorBase<TY> &y){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(x.size() != y.size());
#endif
	TXC sum(0);
	for(size_t i = 0; i < xc.size(); ++i){
		sum += std::conj(xc[i])*y[i];
	}
	return sum;
}

#endif // USE_COMPLEX_MATRICES

//// Scale

template <class T>
MatrixOpStatus Scale(TMatrixBase<T> &A, const T &scale){
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			A(i,j) *= scale;
		}
	}
	return OK;
}
template <class T>
MatrixOpStatus Scale(TVectorBase<T> &x, const T &scale){
	for(size_t i = 0; i < x.size(); ++i){
		x[i] *= scale;
	}
	return OK;
}

//// Add

template <class TB, class TA>
MatrixOpStatus Add(const TMatrixBase<TB> &B, TMatrixBase<TA> &APlusB, const TA &scaleB = TA(1)){
#ifdef USE_MATRIX_OPS_CHECKS
	if(B.Rows() != APlusB.Rows() && B.Cols() != APlusB.Cols()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(B.Rows() != APlusB.Rows());
	assert(B.Cols() != APlusB.Cols());
#endif
	for(size_t j = 0; j < APlusB.Cols(); ++j){
		for(size_t i = 0; i < APlusB.Rows(); ++i){
			APlusB(i,j) += scaleB * B(i,j);
		}
	}
	return OK;
}

template <class TB, class TA>
MatrixOpStatus Add(const TVectorBase<TB> &B, TVectorBase<TA> &APlusB, const TA &scaleB = TA(1)){
#ifdef USE_MATRIX_OPS_CHECKS
	if(B.size() != APlusB.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(B.size() != APlusB.size());
#endif
	for(size_t i = 0; i < APlusB.size(); ++i){
		APlusB[i] += scaleB * B[i];
	}
	return OK;
}

//// Rank1Update
template <class TA, class TX, class TYt>
MatrixOpStatus Rank1Update(TMatrixBase<TA> &A, const TVectorBase<TX> &X, const TVectorBase<TYt> &Yt, const TA &scaleXY = TA(1)){
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != X.size() || A.Cols() != Yt.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Rows() != X.size());
	assert(A.Cols() != Yt.size());
#endif
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			A(i,j) += scaleXY * X[i]*Yt[j];
		}
	}
	return OK;
}

#ifdef USE_COMPLEX_MATRICES

template <class TA, class TX>
MatrixOpStatus Rank1Update(TMatrixBase<TA> &A, const TVectorBase<TX> &X, const TA& scaleXX = TA(1)){
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != X.size() || A.Cols() != X.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Rows() != X.size());
	assert(A.Cols() != X.size());
#endif
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			A(i,j) += scaleXX * X[i]*std::conj(X[j]);
		}
	}
	return OK;
}

#endif // USE_COMPLEX_MATRICES

//// Rank2Update

#ifdef USE_COMPLEX_MATRICES

template <class TA, class TX, class TY>
MatrixOpStatus Rank2Update(TMatrixBase<TA> &A, const TVectorBase<TX> &X, const TVectorBase<TY> &Y, const TA &scaleXY = TA(1)){
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != A.Cols() || A.Rows() != X.size() || A.Cols() != Y.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Rows() != A.Cols());
	assert(A.Rows() != X.size());
	assert(A.Cols() != Y.size());
#endif
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			A(i,j) += scaleXY * (X[i]*std::conj(Y[j]) + Y[i]*std::conj(X[i]));
		}
	}
	return OK;
}

#endif // USE_COMPLEX_MATRICES

//// Mult

template <class TA, class TX, class TY>
MatrixOpStatus Mult(
	const TMatrixBase<TA> &A, const TVectorBase<TX> &X, TVectorBase<TY> &Y,
	const TY &scale_AX = TY(1),
	const TY &scale_Y = TY(0)
){
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Cols() != X.size() || A.Rows() == Y.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
#endif
	for(size_t i = 0; i < A.Rows(); ++i){
		TY sum(0);
		for(size_t j = 0; j < A.Cols(); ++j){
			sum += A(i,j)*X[j];
		}
		Y[i] = scale_AX * sum + scale_Y * Y[i];
	}
	return OK;
}
template <class TA, class TB, class TC>
MatrixOpStatus Mult(
	const TMatrixBase<TA> &A, const TMatrixBase<TB> &B, TMatrixBase<TC> &C,
	const TC &scale_AB = TC(1),
	const TC &scale_C = TC(0)
){
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Cols() != B.Rows() || A.Rows() == C.Rows() || B.Cols() == C.Cols()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
#endif
	for(size_t i = 0; i < A.Rows(); ++i){
		for(size_t j = 0; j < B.Cols(); ++j){
			TC sum(0);
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
inline MatrixOpStatus Mult<
	TMatrix<double>,
	TMatrix<double>,
	TMatrix<double>
>(
	const TMatrix<double> &A, const TMatrix<double> &B, TMatrix<double> &C,
	const double &scale_AB = 1.0,
	const double &scale_C = 0.0
){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
#endif
	TBLAS::TBLAS_NAME(gemm,GEMM)<double,double,double,double,double>(
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
	const double &scale_AB = 1.0,
	const double &scale_C = 0.0
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
