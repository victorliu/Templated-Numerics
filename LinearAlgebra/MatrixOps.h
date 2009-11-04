#ifndef _MATRIX_OPS_H_
#define _MATRIX_OPS_H_

// Preprocessor flags:
//   USE_MATRIX_OPS_CHECKS  - If not turned on, always returns OK
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

template <class GeneralMatrixTypeSrc, class GeneralMatrixTypeDst>
MatrixOpStatus Copy(const typename GeneralMatrixTypeSrc::matrix_type &src, typename GeneralMatrixTypeDst::matrix_type &dst){
#ifdef USE_MATRIX_OPS_CHECKS
	if(src.Rows() != dst.Rows() && src.Cols() != dst.Cols()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(src.Rows() != dst.Rows());
	assert(src.Cols() != dst.Cols());
#endif
	typedef typename GeneralMatrixTypeDst::value_type dst_t;
	for(size_t j = 0; j < dst.Cols(); ++j){
		for(size_t i = 0; i < dst.Rows(); ++i){
			dst(i,j) = dst_t(src(i,j));
		}
	}
	return OK;
}
template <class GeneralVectorTypeSrc, class GeneralVectorTypeDst>
MatrixOpStatus Copy(const typename GeneralVectorTypeSrc::vector_type &src, typename GeneralVectorTypeDst::vector_type &dst){
#ifdef USE_MATRIX_OPS_CHECKS
	if(src.size() != dst.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(src.size() != dst.size());
#endif
	typedef typename GeneralVectorTypeDst::value_type dst_t;
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] = dst_t(src[i]);
	}
	return OK;
}

//// Fill

template <class GeneralMatrixType>
MatrixOpStatus Fill(typename GeneralMatrixType::matrix_type &dst, const typename GeneralMatrixType::value_type &value){
	for(size_t j = 0; j < dst.Cols(); ++j){
		for(size_t i = 0; i < dst.Rows(); ++i){
			dst(i,j) = value;
		}
	}
	return OK;
}
template <class GeneralVectorType>
MatrixOpStatus Fill(typename GeneralVectorType::vector_type &dst, const typename GeneralVectorType::value_type &value){
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] = value;
	}
	return OK;
}

//// Dot

template <class GeneralVectorType1, class GeneralVectorType2>
typename GeneralVectorType1::value_type Dot(const typename GeneralVectorType1::vector_type &x, typename GeneralVectorType2::vector_type &y){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(x.size() != y.size());
#endif
	typename GeneralVectorType1::value_type sum(0);
	for(size_t i = 0; i < dst.Rows(); ++i){
		sum += x[i]*y[i];
	}
	return sum;
}
template <class GeneralVectorType1, class GeneralVectorType2>
typename GeneralVectorType2::value_type Dot(const typename GeneralVectorType1::vector_type &x, typename GeneralVectorType2::vector_type &y){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(x.size() != y.size());
#endif
	typename GeneralVectorType2::value_type sum(0);
	for(size_t i = 0; i < dst.Rows(); ++i){
		sum += x[i]*y[i];
	}
	return sum;
}

//// ConjugateDot

template <class GeneralVectorType1, class GeneralVectorType2>
typename GeneralVectorType1::value_type ConjugateDot(const typename GeneralVectorType1::vector_type &xc, typename GeneralVectorType2::vector_type &y){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(x.size() != y.size());
#endif
	typename GeneralVectorType1::value_type sum(0);
	for(size_t i = 0; i < dst.Rows(); ++i){
		sum += std::conj(x[i])*y[i];
	}
	return sum;
}
template <class GeneralVectorType1, class GeneralVectorType2>
typename GeneralVectorType2::value_type ConjugateDot(const typename GeneralVectorType1::vector_type &xc, typename GeneralVectorType2::vector_type &y){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(x.size() != y.size());
#endif
	typename GeneralVectorType2::value_type sum(0);
	for(size_t i = 0; i < dst.Rows(); ++i){
		sum += std::conj(x[i])*y[i];
	}
	return sum;
}

//// Scale

template <class GeneralMatrixType>
MatrixOpStatus Scale(typename GeneralMatrixType::matrix_type &A, const typename GeneralMatrixType::value_type &scale){
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			A(i,j) *= scale;
		}
	}
	return OK;
}
template <class GeneralVectorType>
MatrixOpStatus Scale(typename GeneralVectorType::vector_type &x, const typename GeneralVectorType::value_type &scale){
	for(size_t i = 0; i < x.size(); ++i){
		x[i] *= scale;
	}
	return OK;
}

//// Add

template <class GeneralMatrixType1, class GeneralMatrixType2>
MatrixOpStatus Add(const typename GeneralMatrixType1::matrix_type &B, typename GeneralMatrixType2::matrix_type &APlusB, const typename GeneralMatrixType2::value_type &scaleB = typename GeneralMatrixType2::value_type(1)){
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

template <class GeneralVectorType1, class GeneralVectorType2>
MatrixOpStatus Add(const typename GeneralVectorType1::vector_type &B, typename GeneralVectorType2::vector_type &APlusB, const typename GeneralVectorType2::value_type &scaleB = typename GeneralVectorType2::value_type(1)){
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
MatrixOpStatus Rank1Update(typename TA::matrix_type &A, const typename TX::vector_type &X, const typename TYt::vector_type &Yt, const TA::value_type &scaleXY = typename TA::value_type(1)){
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

template <class TA, class TX, class TYt>
MatrixOpStatus Rank1Update(typename TA::matrix_type &A, const typename TX::vector_type &X, const TA::value_type &scaleXX = typename TA::value_type(1)){
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

//// Rank2Update

template <class TA, class TX, class TY>
MatrixOpStatus Rank2Update(typename TA::matrix_type &A, const typename TX::vector_type &X, const typename TY::vector_type &Y, const TA::value_type &scaleXY = typename TA::value_type(1)){
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
			A(i,j) += scaleXY * (X[i]*std::conj(Y[j]) + Y[i]*std::conj(X[i]);
		}
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
inline MatrixOpStatus Mult<
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
