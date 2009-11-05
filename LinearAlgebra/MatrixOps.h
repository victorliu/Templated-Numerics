#ifndef _MATRIX_OPS_H_
#define _MATRIX_OPS_H_

// Preprocessor flags:
//   USE_MATRIX_OPS_CHECKS  - If not turned on, always returns OK
//   USE_MATRIX_ASSERTS - Turn on assertion checking
//   USE_MATRIX_OPS_TBLAS   - Use optimized BLAS routines when possible
//   USE_COMPLEX_MATRICES
//   USE_ADVANCED_MATRIX_OPS

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
//   Swap(x, y)
//   LargestElementIndex(X)
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
//   LUDecomposition(A, P) - A -> P*unit_lower(A)*upper(A)
//   Solve(A, B, X) - A*X == B
//   SolveLeastSquares(A, B, X) - A*X == B
//   Eigenvalues(A, L)
//   Eigensystem(A, L, X) - A*X == X*diag(L)
//   SingularValueDecomposition(A, U, S, V) - A == U*S*V
//   CholeskyDecomposition(A, L) - A == L*L'

//// Copy

template <class TSRC, class TDST>
MatrixOpStatus Copy(const TMatrixBase<TSRC> &src, TMatrixBase<TDST> &dst){
#ifdef USE_MATRIX_OPS_CHECKS
	if(src.Rows() != dst.Rows() && src.Cols() != dst.Cols()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
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
	assert(src.size() == dst.size());
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

//// Swap
template <class T>
MatrixOpStatus Swap(TVectorBase<T> &x, TVectorBase<T> &y){
#ifdef USE_MATRIX_OPS_CHECKS
	if(x.size() != y.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(x.size() == y.size());
#endif
	for(size_t i = 0; i < x.size(); ++i){
		std::swap(x[i], y[i]);
	}
	return OK;
}
template <class ViewType>
MatrixOpStatus Swap(ViewType x, ViewType y){
#ifdef USE_MATRIX_OPS_CHECKS
	if(x.size() != y.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(x.size() == y.size());
#endif
	for(size_t i = 0; i < x.size(); ++i){
		std::swap(x[i], y[i]);
	}
	return OK;
}

//// LargestElementIndex

template <class T, template <class TT> class TV>
size_t LargestElementIndex(const TVectorBase<TV<T> > &x){
	size_t ret = 0;
	T mag = std::abs(x[0]);
	for(size_t i = 1; i < x.size(); ++i){
		T new_mag = std::abs(x[i]);
		if(new_mag > mag){
			ret = i;
			mag = new_mag;
		}
	}
	return ret;
}
template <class T>
size_t LargestElementIndex(const TVectorBase<T> &x){
	size_t ret = 0;
	T mag = std::abs(x[0]);
	for(size_t i = 1; i < x.size(); ++i){
		T new_mag = std::abs(x[i]);
		if(new_mag > mag){
			ret = i;
			mag = new_mag;
		}
	}
	return ret;
}

//// Dot

template <class TX, class TY>
TX Dot(const TVectorBase<TX> &x, const TVectorBase<TY> &y){
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(x.size() == y.size());
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
	assert(xc.size() == y.size());
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
template <class ViewType>
MatrixOpStatus Scale(ViewType A, const typename ViewType::value_type &scale){
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			A(i,j) *= scale;
		}
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
	assert(B.Rows() == APlusB.Rows());
	assert(B.Cols() == APlusB.Cols());
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
	assert(B.size() == APlusB.size());
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
	assert(A.Rows() == X.size());
	assert(A.Cols() == Yt.size());
#endif
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			A(i,j) += scaleXY * X[i]*Yt[j];
		}
	}
	return OK;
}
template <class ViewA, class ViewX, class ViewYt>
MatrixOpStatus Rank1Update(ViewA A, ViewX X, ViewYt Yt, const typename ViewA::value_type &scaleXY = typename ViewA::value_type(1)){
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != X.size() || A.Cols() != Yt.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Rows() == X.size());
	assert(A.Cols() == Yt.size());
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
	assert(A.Rows() == X.size());
	assert(A.Cols() == X.size());
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
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == X.size());
	assert(A.Cols() == Y.size());
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






















#ifdef USE_ADVANCED_MATRIX_OPS

//// LUDecomposition(A, P, L, U) - A == P*L*U
template <class TA, class TP>
MatrixOpStatus LUDecomposition(TMatrixBase<TA> &A, TVectorBase<TP> &Pivots){
	const size_t min_dim = ((A.Rows() < A.Cols()) ? A.Rows() : A.Cols());
#ifdef USE_MATRIX_OPS_CHECKS
	if(min_dim != Pivots.size()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(min_dim == Pivots.size());
#endif
	size_t info = 0;
	for(size_t j = 0; j < min_dim; ++j){
		size_t jp = j + LargestElementIndex(SubVector(GetColumn(A, j), j));
		Pivots[j] = jp;
		if(TA(0) != A(jp,j)){
			if(jp != j){
				Swap(GetRow(A, j), GetRow(A, jp));
			}
			if(j < A.Rows()){
				Scale(SubVector(GetColumn(A,j),j+1,A.Rows()-j-1), TA(1)/A(j,j)); // possible overflow when inverting A(j,j)
			}
		}else{
			info = j;
		}
		if(j < min_dim){
			Rank1Update(SubMatrix(A, j+1,j+1, A.Rows()-j-1, A.Cols()-j-1), SubVector(GetColumn(A,j), j+1), SubVector(GetRow(A,j),j+1), TA(-1));
		}
	}
	if(0 != info){ return SINGULAR_MATRIX; }
	else{ return OK; }
}


//// Solve(A, B, X) - A*X == B
template <class TA, class TB, class TX>
MatrixOpStatus Solve(const TA &A, const TMatrixBase<TB> &B, TMatrixBase<TX> &X){
	TA Acopy(A);
	Copy(B, X);
	TVector<size_t> Pivots(A.Rows());
	LUDecomposition(Acopy, Pivots);
	// Apply pivots
	for(size_t i = 0; i < Pivots.size(); ++i){
		if(Pivots[i] != i){
			Swap(GetRow(X,i), GetRow(X,Pivots[i]));
		}
	}
	// Solve lower unit
	for(size_t j = 0; j < X.Cols(); ++j){
		for(size_t k = 0; k < X.Rows(); ++k){
			if(typename TMatrixBase<TX>::value_type(0) != X(k,j)){
				for(size_t i = k+1; i < X.Rows(); ++i){
					X(i,j) -= X(k,j)*Acopy(i,k);
				}
			}
		}
	}
	// Solver upper non unit
	for(size_t j = 0; j < X.Cols(); ++j){
		for(size_t k = X.Rows()-1; (signed)k >= 0; --k){
			if(typename TMatrixBase<TX>::value_type(0) != X(k,j)){
				X(k,j) /= Acopy(k,k);
				for(size_t i = 0; i < k; ++i){
					X(i,j) -= X(k,j)*Acopy(i,k);
				}
			}
		}
	}
}
template <class TA, class TB, class TX>
MatrixOpStatus Solve(const TA &A, const TVectorBase<TB> &B, TVectorBase<TX> &X){
	TA Acopy(A);
	Copy(B, X);
	TVector<size_t> Pivots(A.Rows());
	LUDecomposition(Acopy, Pivots);
	// Apply pivots
	for(size_t i = 0; i < Pivots.size(); ++i){
		if(Pivots[i] != i){
			Swap(GetRow(X,i), GetRow(X,Pivots[i]));
		}
	}
	// Solve lower unit
	for(size_t k = 0; k < X.size(); ++k){
		if(typename TVectorBase<TX>::value_type(0) != X[k]){
			for(size_t i = k+1; i < X.size(); ++i){
				X[i] -= X[k]*Acopy(i,k);
			}
		}
	}
	// Solver upper non unit
	for(size_t k = X.size()-1; (signed)k >= 0; --k){
		if(typename TVectorBase<TX>::value_type(0) != X[k]){
			X[k] /= Acopy(k,k);
			for(size_t i = 0; i < k; ++i){
				X[i] -= X[k]*Acopy(i,k);
			}
		}
	}
}

#endif // USE_ADVANCED_MATRIX_OPS

}; // namespace MatrixOps

#endif // _MATRIX_OPS_H_
