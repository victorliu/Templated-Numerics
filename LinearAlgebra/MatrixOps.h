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

#include "TMatrix.h"
#include "MatrixViews.h"
#include "TDiagonalMatrix.h"

#ifdef USE_MATRIX_ASSERTS
# include <cassert>
#endif

#ifdef USE_COMPLEX_MATRICES
# include <complex>
#endif

#ifdef USE_MATRIX_OPS_TBLAS
# include "TBLAS/tblas.hpp"
#endif

#ifndef MATRIX_OP_STATUS_DEFINED
#define MATRIX_OP_STATUS_DEFINED
enum MatrixOpStatus{
	OK,
	DIMENSION_MISMATCH,
	SINGULAR_MATRIX,
	UNKNOWN_ERROR
};
#endif

namespace MatrixOps{

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
//   Mult(D, A); - D is diagonal
//   Mult(A, D); - D is diagonal
//   FrobeniusNorm(A)
//
// Sophisticated operations
//   LUDecomposition(A, P) - A -> P*unit_lower(A)*upper(A)
//   SolveDestructive(A, X) - X holds B, A gets overwritten with LU, X overwritten with solution
//   Solve(A, B, X) - A*X == B
//   Invert(A)
//   Invert(A, Ainv)
//   SolveLeastSquares(A, B, X) - A*X == B
//   Eigenvalues(A, L)
//   Eigensystem(A, L, X) - A*X == X*diag(L)
//   SingularValueDecomposition(A, U, S, V) - A == U*S*V
//   CholeskyDecomposition(A, L) - A == L*L'
//   UnitaryProcrustes(A) - Replaces A with the nearest (Frobenius norm) unitary matrix, A'*A = I
//   GeneralizedProcrustes(A, B, C) - Replaces A with nearest (Frobenius norm) matrix such that A'*B*A = C


//// Copy

template <class Tsrc, class Tdst>
	typename IsReadableMatrix<typename Tsrc::readable_matrix,
	typename IsWritableMatrixView<typename Tdst::writable_matrix,
MatrixOpStatus
	>::type>::type
Copy(const Tsrc &src, const Tdst &dst){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
	typedef typename Tdst::value_type dest_type;
#ifdef USE_MATRIX_OPS_CHECKS
	if(src.Rows() != dst.Rows() && src.Cols() != dst.Cols()){ return DIMENSION_MISMATCH; }
#endif
	
	for(size_t j = 0; j < dst.Cols(); ++j){
		for(size_t i = 0; i < dst.Rows(); ++i){
			dst(i,j) = dest_type(src(i,j));
		}
	}
	return OK;
}
template <class Tsrc, class Tdst>
	typename IsReadableMatrix<typename Tsrc::readable_matrix,
	typename IsWritableMatrix<typename Tdst::writable_matrix,
MatrixOpStatus
	>::type>::type
Copy(const Tsrc &src, Tdst &dst){
	return Copy(src, TrivialWritableMatrixView<Tdst>(dst));
}


template <class Tsrc, class Tdst>
	typename IsReadableVector<typename Tsrc::readable_vector,
	typename IsWritableVectorView<typename Tdst::writable_vector,
MatrixOpStatus
	>::type>::type
Copy(const Tsrc &src, const Tdst &dst){
	assert(src.size() == dst.size());
	typedef typename Tdst::value_type dest_type;
#ifdef USE_MATRIX_OPS_CHECKS
	if(src.size() != dst.size()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] = dest_type(src[i]);
	}
	return OK;
}
template <class Tsrc, class Tdst>
	typename IsReadableVector<typename Tsrc::readable_vector,
	typename IsWritableVector<typename Tdst::writable_vector,
MatrixOpStatus
	>::type>::type
Copy(const Tsrc &src, Tdst &dst){
	return Copy(src, TrivialWritableVectorView<Tdst>(dst));
}

//// Fill

template <class Tdst>
	typename IsWritableMatrixView<typename Tdst::writable_matrix,
MatrixOpStatus
	>::type
Fill(const Tdst &dst, const typename Tdst::value_type &value){
	for(size_t j = 0; j < dst.Cols(); ++j){
		for(size_t i = 0; i < dst.Rows(); ++i){
			dst(i,j) = value;
		}
	}
	return OK;
}
template <class Tdst>
	typename IsWritableMatrix<typename Tdst::writable_matrix,
MatrixOpStatus
	>::type
Fill(Tdst &dst, const typename Tdst::value_type &value){
	return Fill(TrivialWritableMatrixView<Tdst>(dst), value);
}

template <class Tdst>
	typename IsWritableVectorView<typename Tdst::writable_vector,
MatrixOpStatus
	>::type
Fill(const Tdst &dst, const typename Tdst::value_type &value){
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] = value;
	}
	return OK;
}
template <class Tdst>
	typename IsWritableVector<typename Tdst::writable_vector,
MatrixOpStatus
	>::type
Fill(Tdst &dst, const typename Tdst::value_type &value){
	return Fill(TrivialWritableVectorView<Tdst>(dst), value);
}


//// Swap

template <class TX, class TY>
	typename IsWritableVectorView<typename TX::writable_vector,
	typename IsWritableVectorView<typename TY::writable_vector,
MatrixOpStatus
	>::type>::type
Swap(const TX &x, const TY &y){
	assert(x.size() == y.size());
#ifdef USE_MATRIX_OPS_CHECKS
	if(x.size() != y.size()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t i = 0; i < x.size(); ++i){
		std::swap(x[i], y[i]);
	}
	return OK;
}
template <class TX, class TY>
	typename IsWritableVector<typename TX::writable_vector,
	typename IsWritableVectorView<typename TY::writable_vector,
MatrixOpStatus
	>::type>::type
Swap(TX &x, const TY &y){
	return Swap(TrivialWritableVectorView<TX>(x), y);
}
template <class TX, class TY>
	typename IsWritableVectorView<typename TX::writable_vector,
	typename IsWritableVector<typename TY::writable_vector,
MatrixOpStatus
	>::type>::type
Swap(const TX &x, TY &y){
	return Swap(x, TrivialWritableVectorView<TY>(y));
}
template <class TX, class TY>
	typename IsWritableVector<typename TX::writable_vector,
	typename IsWritableVector<typename TY::writable_vector,
MatrixOpStatus
	>::type>::type
Swap(const TX &x, TY &y){
	return Swap(TrivialWritableVectorView<TX>(x), TrivialWritableVectorView<TY>(y));
}



//// LargestElementIndex

template <class T>
	typename IsReadableVector<typename T::readable_vector,
size_t
	>::type
LargestElementIndex(const T &x){
	typedef typename T::value_type value_Type;
	size_t ret = 0;
	for(size_t i = 1; i < x.size(); ++i){
		if(std::abs(x[i]) > std::abs(x[ret])){
			ret = i;
		}
	}
	return ret;
}

//// Dot

template <class TX, class TY>
	typename IsReadableVector<typename TX::readable_vector,
	typename IsReadableVector<typename TY::readable_vector,
typename TX::value_type
	>::type>::type
Dot(const TX &x, const TY &y){
	assert(x.size() == y.size());
	typename TX::value_type sum(0);
	for(size_t i = 0; i < x.size(); ++i){
		sum += x[i]*y[i];
	}
	return sum;
}

//// ConjugateDot

#ifdef USE_COMPLEX_MATRICES

template <class TXC, class TY>
	typename IsReadableVector<typename TXC::readable_vector,
	typename IsReadableVector<typename TY::readable_vector,
typename TXC::value_type
	>::type>::type
ConjugateDot(const TXC &xc, const TY &y){
	assert(xc.size() == y.size());
	typename TXC::value_type sum(0);
	for(size_t i = 0; i < xc.size(); ++i){
		sum += std::conj(xc[i])*y[i];
	}
	return sum;
}

#endif // USE_COMPLEX_MATRICES

//// Scale

template <class TA>
	typename IsWritableMatrixView<typename TA::writable_matrix,
MatrixOpStatus
	>::type
Scale(const TA &A, const typename TA::value_type &scale){
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			A(i,j) *= scale;
		}
	}
	return OK;
}
template <class TA>
	typename IsWritableMatrix<typename TA::writable_matrix,
MatrixOpStatus
	>::type
Scale(TA &A, const typename TA::value_type &scale){
	return Scale(TrivialWritableMatrixView<TA>(A), scale);
}

template <class TX>
	typename IsWritableVectorView<typename TX::writable_vector,
MatrixOpStatus
	>::type
Scale(const TX &x, const typename TX::value_type &scale){
	for(size_t i = 0; i < x.size(); ++i){
		x[i] *= scale;
	}
	return OK;
}
template <class TX>
	typename IsWritableVector<typename TX::writable_vector,
MatrixOpStatus
	>::type
Scale(TX &x, const typename TX::value_type &scale){
	return Scale(TrivialWritableVectorView<TX>(x), scale);
}

//// Add

template <class Tsrc, class Tdst>
	typename IsReadableMatrix<typename Tsrc::readable_matrix,
	typename IsWritableMatrixView<typename Tdst::writable_matrix,
MatrixOpStatus
	>::type>::type
Add(const Tsrc &src, const Tdst &dst, const typename Tdst::value_type &scale_src = typename Tdst::value_type(1)){
	assert(src.Rows() == dst.Rows());
	assert(src.Cols() == dst.Cols());
#ifdef USE_MATRIX_OPS_CHECKS
	if(src.Rows() != dst.Rows() && src.Cols() != dst.Cols()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t j = 0; j < dst.Cols(); ++j){
		for(size_t i = 0; i < dst.Rows(); ++i){
			dst(i,j) += scale_src * src(i,j);
		}
	}
	return OK;
}
template <class Tsrc, class Tdst>
	typename IsReadableMatrix<typename Tsrc::readable_matrix,
	typename IsWritableMatrix<typename Tdst::writable_matrix,
MatrixOpStatus
	>::type>::type
Add(const Tsrc &src, Tdst &dst, const typename Tdst::value_type &scale_src = typename Tdst::value_type(1)){
	return Add(src, TrivialWritableMatrixView<Tdst>(dst), scale_src);
}

template <class Tsrc, class Tdst>
	typename IsReadableVector<typename Tsrc::readable_vector,
	typename IsWritableVectorView<typename Tdst::writable_vector,
MatrixOpStatus
	>::type>::type
Add(const Tsrc &src, const Tdst &dst, const typename Tdst::value_type &scale_src = typename Tdst::value_type(1)){
	assert(src.size() == dst.size());
#ifdef USE_MATRIX_OPS_CHECKS
	if(src.size() != dst.size()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t i = 0; i < dst.size(); ++i){
		dst[i] += scale_src * src[i];
	}
	return OK;
}
template <class Tsrc, class Tdst>
	typename IsReadableVector<typename Tsrc::readable_vector,
	typename IsWritableVector<typename Tdst::writable_vector,
MatrixOpStatus
	>::type>::type
Add(const Tsrc &src, Tdst &dst, const typename Tdst::value_type &scale_src = typename Tdst::value_type(1)){
	return Add(src, TrivialWritableVectorView<Tdst>(dst), scale_src);
}


//// Rank1Update

template <class TA, class TX, class TYt>
	typename IsWritableMatrixView<typename TA::writable_matrix,
	typename IsReadableVector<typename TX::readable_vector,
	typename IsReadableVector<typename TYt::readable_vector,
MatrixOpStatus
	>::type>::type>::type
Rank1Update(const TA &A, const TX &X, const TYt &Yt, const typename TA::value_type &scale_XY = typename TA::value_type(1)){
	assert(A.Rows() == X.size());
	assert(A.Cols() == Yt.size());
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != X.size() || A.Cols() != Yt.size()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			A(i,j) += scale_XY * X[i]*Yt[j];
		}
	}
	return OK;
}
template <class TA, class TX, class TYt>
	typename IsWritableMatrix<typename TA::writable_matrix,
	typename IsReadableVector<typename TX::readable_vector,
	typename IsReadableVector<typename TYt::readable_vector,
MatrixOpStatus
	>::type>::type>::type
Rank1Update(TA &A, const TX &X, const TYt &Yt, const typename TA::value_type &scale_XY = typename TA::value_type(1)){
	return Rank1Update(TrivialWritableMatrixView<TA>(A), X, Yt, scale_XY);
}
/*
#ifdef USE_COMPLEX_MATRICES

template <class TA, class TX>
MatrixOpStatus Rank1Update(const WritableMatrixView<TA> &A, const ReadableVector<TX> &X, const TA& scaleXX = TA(1)){
	assert(A.Rows() == X.size());
	assert(A.Cols() == X.size());
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != X.size() || A.Cols() != X.size()){ return DIMENSION_MISMATCH; }
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
MatrixOpStatus Rank2Update(const WritableMatrixView<TA> &A, const ReadableVector<TX> &X, const ReadableVector<TY> &Y, const TA &scaleXY = TA(1)){
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

*/

//// Mult matrix vector

template <class TA, class TX, class TY>
	typename IsReadableMatrix<typename TA::readable_matrix,
	typename IsReadableVector<typename TX::readable_vector,
	typename IsWritableVectorView<typename TY::writable_vector,
MatrixOpStatus
	>::type>::type>::type
Mult(const TA &A, const TX &X, const TY &Y, const typename TY::value_type &scale_AX = typename TY::value_type(1), const typename TY::value_type &scale_Y = typename TY::value_type(0)){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Cols() != X.size() || A.Rows() == Y.size()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t i = 0; i < A.Rows(); ++i){
		typename TY::value_type sum(0);
		for(size_t j = 0; j < A.Cols(); ++j){
			sum += A(i,j)*X[j];
		}
		Y[i] = scale_AX * sum + scale_Y * Y[i];
	}
	return OK;
}
template <class TA, class TX, class TY>
	typename IsReadableMatrix<typename TA::readable_matrix,
	typename IsReadableVector<typename TX::readable_vector,
	typename IsWritableVector<typename TY::writable_vector,
MatrixOpStatus
	>::type>::type>::type
Mult(const TA &A, const TX &X, TY &Y, const typename TY::value_type &scale_AX = typename TY::value_type(1), const typename TY::value_type &scale_Y = typename TY::value_type(0)){
	return Mult(A, X, TrivialWritableVectorView<TY>(Y), scale_AX, scale_Y);
}


//// Mult matrix matrix
template <class TA, class TB, class TC>
	typename IsReadableMatrix<typename TA::readable_matrix,
	typename IsReadableMatrix<typename TB::readable_matrix,
	typename IsWritableMatrixView<typename TC::writable_matrix,
MatrixOpStatus
	>::type>::type>::type
Mult(const TA &A, const TB &B, const TC &C, const typename TC::value_type &scale_AB = typename TC::value_type(1), const typename TC::value_type &scale_C = typename TC::value_type(0)){
	assert(A.Cols() == B.Rows());
	assert(A.Rows() == C.Rows());
	assert(B.Cols() == C.Cols());
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Cols() != B.Rows() || A.Rows() == C.Rows() || B.Cols() == C.Cols()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t i = 0; i < A.Rows(); ++i){
		for(size_t j = 0; j < B.Cols(); ++j){
			typename TC::value_type sum(0);
			for(size_t k = 0; k < A.Cols(); ++k){
				sum += A(i,k)*B(k,j);
			}
			C(i,j) = scale_AB * sum + scale_C * C(i,j);
		}
	}
	return OK;
}
template <class TA, class TB, class TC>
	typename IsReadableMatrix<typename TA::readable_matrix,
	typename IsReadableMatrix<typename TB::readable_matrix,
	typename IsWritableMatrix<typename TC::writable_matrix,
MatrixOpStatus
	>::type>::type>::type
Mult(const TA &A, const TB &B, TC &C, const typename TC::value_type &scale_AB = typename TC::value_type(1), const typename TC::value_type &scale_C = typename TC::value_type(0)){
	return Mult(A, B, TrivialWritableMatrixView<TC>(C), scale_AB, scale_C);
}




//// Mult diagonal

template <class TD, class TA>
	typename IsWritableMatrixView<typename TA::writable_matrix,
MatrixOpStatus
	>::type
Mult(const TDiagonalMatrix<TD> &D, const TA &A){
	assert(D.size() == A.Rows());
#ifdef USE_MATRIX_OPS_CHECKS
	if(D.size() != A.Rows()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t i = 0; i < A.Rows(); ++i){
		Scale(GetRow(A,i), D[i]);
	}
	return OK;
}
template <class TD, class TA>
	typename IsWritableMatrix<typename TA::writable_matrix,
MatrixOpStatus
	>::type
Mult(const TDiagonalMatrix<TD> &D, TA &A){
	return Mult(D, TrivialWritableMatrixView<TA>(A));
}

template <class TA, class TD>
	typename IsWritableMatrixView<typename TA::writable_matrix,
MatrixOpStatus
	>::type
Mult(const TA &A, const TDiagonalMatrix<TD> &D){
	assert(D.size() == A.Cols());
#ifdef USE_MATRIX_OPS_CHECKS
	if(D.size() != A.Cols()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t j = 0; j < A.Cols(); ++j){
		Scale(GetColumn(A,j), D[j]);
	}
	return OK;
}
template <class TA, class TD>
	typename IsWritableMatrix<typename TA::writable_matrix,
MatrixOpStatus
	>::type
Mult(TA &A, const TDiagonalMatrix<TD> &D){
	return Mult(TrivialWritableMatrixView<TA>(A), D);
}


template <class TD, class TX>
	typename IsWritableVectorView<typename TX::writable_vector,
MatrixOpStatus
	>::type
Mult(const TDiagonalMatrix<TD> &D, const TX &X){
	assert(D.size() == X.size());
#ifdef USE_MATRIX_OPS_CHECKS
	if(D.size() != X.size()){ return DIMENSION_MISMATCH; }
#endif
	for(size_t i = 0; i < X.size(); ++i){
		X[i] *= D[i];
	}
	return OK;
}
template <class TD, class TX>
	typename IsWritableVector<typename TX::writable_vector,
MatrixOpStatus
	>::type
Mult(const TDiagonalMatrix<TD> &D, TX &X){
	return Mult(D, TrivialWritableMatrixView<TX>(X));
}



#ifdef USE_COMPLEX_MATRICES
template <class TA>
	typename IsReadableMatrix<typename TA::readable_matrix,
typename TA::value_type
	>::type
FrobeniusNorm(const TA &A){
	typename TA::value_type sum(0);
	for(size_t j = 0; j < A.Cols(); ++j){
		for(size_t i = 0; i < A.Rows(); ++i){
			typename TA::value_type t(std::abs(A(i,j)));
			sum += t*t;
		}
	}
	return sum;
}
#endif // USE_COMPLEX_MATRICES




#ifdef USING_TCCS_MATRIX

#include "TCCSMatrix.h"

template <class TA, class TX, class TY, class TAlloc>
	typename IsWritableVector<typename TY::writable_vector,
MatrixOpStatus
	>::type
Mult(
	const TCCSMatrix<TA,TAlloc> &A, const TX &X, TY &Y,
	const typename TY::value_type &scale_AX = typename TY::value_type(1),
	const typename TY::value_type &scale_Y = typename TY::value_type(0)
){
	assert(A.Cols() == X.size());
	assert(A.Rows() == Y.size());
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Cols() != X.size() || A.Rows() == Y.size()){ return DIMENSION_MISMATCH; }
#endif
	Scale(Y, scale_Y);
	if(A.flags & TCCSMatrix<TA,TAlloc>::SYMMETRIC){
		for(int j = 0; j < (int)A.Cols(); j++){
			for(int ip = A.colptr[j]; ip < A.colptr[j+1]; ip++){
				int i = A.rowind[ip];
				Y[i] += scale_AX*X[j]*A.values[ip];
				if(i != j){ Y[j] += scale_AX*X[i]*A.values[ip]; }
			}
		}
	}
#ifdef USE_COMPLEX_MATRICES
	else if(A.flags & TCCSMatrix<TA,TAlloc>::HERMITIAN){
		for(int j=0; j < (int)A.Cols(); j++){
			for(int ip = A.colptr[j]; ip < A.colptr[j+1]; ip++){
				int i = A.rowind[ip];
				Y[i] += scale_AX*X[j]*A.values[ip];
				if(i != j){ Y[j] += scale_AX*X[i]*std::conj(A.values[ip]); }
			}
		}
	}
#endif	
	else{
		for(int j = 0; j < (int)A.Cols(); j++){
			for(int ip = A.colptr[j]; ip < A.colptr[j+1]; ip++){
				int i = A.rowind[ip];
				Y[i] += scale_AX*X[j]*A.values[ip];
			}
		}
	}
	return OK;
}

template <class TA, class TX, class TY, class TAlloc>
	typename IsWritableVector<typename TY::writable_vector,
MatrixOpStatus
	>::type
Mult(
	const ConjugateTransposeView<TrivialReadableMatrixView<TCCSMatrix<TA,TAlloc> > > &A, const TX &X, TY &Y,
	const typename TY::value_type &scale_AX = typename TY::value_type(1),
	const typename TY::value_type &scale_Y = typename TY::value_type(0)
){
	assert(A.Rows() == X.size());
	assert(A.Cols() == Y.size());
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != X.size() || A.Cols() == Y.size()){ return DIMENSION_MISMATCH; }
#endif
	Scale(Y, scale_Y);
	if(A.flags & TCCSMatrix<TA,TAlloc>::SYMMETRIC){
		for(int j = 0; j < (int)A.Cols(); j++){
			for(int ip = A.colptr[j]; ip < A.colptr[j+1]; ip++){
				int i = A.rowind[ip];
				Y[j] += scale_AX*X[i]*std::conj(A.values[ip]);
				if(i != j){ Y[i] += scale_AX*X[j]*std::conj(A.values[ip]); }
			}
		}
	}
#ifdef USE_COMPLEX_MATRICES
	else if(A.flags & TCCSMatrix<TA,TAlloc>::HERMITIAN){
		for(int j=0; j < (int)A.Cols(); j++){
			for(int ip = A.colptr[j]; ip < A.colptr[j+1]; ip++){
				int i = A.rowind[ip];
				Y[j] += scale_AX*X[i]*A.values[ip];
				if(i != j){ Y[i] += scale_AX*X[j]*std::conj(A.values[ip]); }
			}
		}
	}
#endif
	else{
		for(int j = 0; j < (int)A.Cols(); j++){
			for(int ip = A.colptr[j]; ip < A.colptr[j+1]; ip++){
				int i = A.rowind[ip];
				Y[j] += scale_AX*X[j]*std::conj(A.values[ip]);
			}
		}
	}
	return OK;
}

#endif // USING_TCCS_MATRIX
















#ifdef USE_ADVANCED_MATRIX_OPS

//// LUDecomposition(A, P, L, U) - A == P*L*U

template <class TA>
	typename IsWritableMatrixView<typename TA::writable_matrix,
MatrixOpStatus
	>::type
LUDecomposition(const TA &A, WritableVector<size_t> &Pivots){
	typedef typename TA::value_type value_type;
	const size_t min_dim = ((A.Rows() < A.Cols()) ? A.Rows() : A.Cols());
	assert(min_dim == Pivots.size());
#ifdef USE_MATRIX_OPS_CHECKS
	if(min_dim != Pivots.size()){ return DIMENSION_MISMATCH; }
#endif
	size_t info = 0;
	for(size_t j = 0; j < min_dim; ++j){
		size_t jp = j + LargestElementIndex(SubVector(GetColumn(A, j), j, A.Rows()-j));
		Pivots[j] = jp;
		if(value_type(0) != A(jp,j)){
			if(jp != j){
				Swap(GetRow(A, j), GetRow(A, jp));
			}
			if(j < A.Rows()){
				Scale(SubVector(GetColumn(A,j),j+1,A.Rows()-j-1), value_type(1)/A(j,j)); // possible overflow when inverting A(j,j)
			}
		}else{
			info = j;
		}
		if(j < min_dim){
			Rank1Update(SubMatrix(A, j+1,j+1, A.Rows()-j-1, A.Cols()-j-1), SubVector(GetColumn(A,j), j+1, A.Rows()-j-1), SubVector(GetRow(A,j),j+1, A.Cols()-j-1), value_type(-1));
		}
	}
	if(0 != info){ return SINGULAR_MATRIX; }
	else{ return OK; }
}
template <class TA>
	typename IsWritableMatrix<typename TA::writable_matrix,
MatrixOpStatus
	>::type
LUDecomposition(TA &A, WritableVector<size_t> &Pivots){
	return LUDecomposition(TrivialWritableMatrixView<TA>(A), Pivots);
}


//// SolveDestructive(A, X) - X holds B, A gets overwritten with LU, X overwritten with solution
template <class TA, class TX>
	typename IsWritableMatrixView<typename TA::writable_matrix,
	typename IsWritableMatrixView<typename TX::writable_matrix,
MatrixOpStatus
	>::type>::type
SolveDestructive(const TA &A, const TX &X){
	typedef typename TX::value_type value_type;
	TVector<size_t> Pivots(A.Rows());
	MatrixOpStatus ret;
	ret = LUDecomposition(A, Pivots);
	// Apply pivots
	for(size_t i = 0; i < Pivots.size(); ++i){
		if(Pivots[i] != i){
			Swap(GetRow(X,i), GetRow(X,Pivots[i]));
		}
	}
	// Solve lower unit
	for(size_t j = 0; j < X.Cols(); ++j){
		for(size_t k = 0; k < X.Rows(); ++k){
			if(value_type(0) != X(k,j)){
				for(size_t i = k+1; i < X.Rows(); ++i){
					X(i,j) -= X(k,j)*A(i,k);
				}
			}
		}
	}
	// Solver upper non unit
	for(size_t j = 0; j < X.Cols(); ++j){
		for(size_t k = X.Rows()-1; (signed)k >= 0; --k){
			if(value_type(0) != X(k,j)){
				X(k,j) /= A(k,k);
				for(size_t i = 0; i < k; ++i){
					X(i,j) -= X(k,j)*A(i,k);
				}
			}
		}
	}
	return ret;
}
template <class TA, class TX>
	typename IsWritableMatrix<typename TA::writable_matrix,
	typename IsWritableMatrixView<typename TX::writable_matrix,
MatrixOpStatus
	>::type>::type
SolveDestructive(TA &A, const TX &X){
	return SolveDestructive(TrivialWritableMatrixView<TA>(A), X);
}
template <class TA, class TX>
	typename IsWritableMatrixView<typename TA::writable_matrix,
	typename IsWritableMatrix<typename TX::writable_matrix,
MatrixOpStatus
	>::type>::type
SolveDestructive(const TA &A, TX &X){
	return SolveDestructive(A, TrivialWritableMatrixView<WritableMatrix<TX> >(X));
}
template <class TA, class TX>
	typename IsWritableMatrix<typename TA::writable_matrix,
	typename IsWritableMatrix<typename TX::writable_matrix,
MatrixOpStatus
	>::type>::type
SolveDestructive(TA &A, TX &X){
	return SolveDestructive(TrivialWritableMatrixView<TA>(A), TrivialWritableMatrixView<TX>(X));
}



template <class TA, class TX>
	typename IsWritableMatrixView<typename TA::writable_matrix,
	typename IsWritableVectorView<typename TX::writable_vector,
MatrixOpStatus
	>::type>::type
SolveDestructive(const TA &A, const TX &X){
	typedef typename TX::value_type value_type;
	TVector<size_t> Pivots(A.Rows());
	MatrixOpStatus ret;
	ret = LUDecomposition(A, Pivots);
	// Apply pivots
	for(size_t i = 0; i < Pivots.size(); ++i){
		if(Pivots[i] != i){
			std::swap(X[i], X[Pivots[i]]);
		}
	}
	// Solve lower unit
	for(size_t k = 0; k < X.size(); ++k){
		if(value_type(0) != X[k]){
			for(size_t i = k+1; i < X.size(); ++i){
				X[i] -= X[k]*A(i,k);
			}
		}
	}
	// Solver upper non unit
		for(size_t k = X.size()-1; (signed)k >= 0; --k){
			if(value_type(0) != X[k]){
				X[k] /= A(k,k);
				for(size_t i = 0; i < k; ++i){
					X[i] -= X[k]*A(i,k);
				}
			}
		}
	return ret;
}
template <class TA, class TX>
	typename IsWritableMatrix<typename TA::writable_matrix,
	typename IsWritableVectorView<typename TX::writable_vector,
MatrixOpStatus
	>::type>::type
SolveDestructive(TA &A, const TX &X){
	return SolveDestructive(TrivialWritableMatrixView<TA>(A), X);
}
template <class TA, class TX>
	typename IsWritableMatrixView<typename TA::writable_matrix,
	typename IsWritableVector<typename TX::writable_vector,
MatrixOpStatus
	>::type>::type
SolveDestructive(const TA &A, TX &X){
	return SolveDestructive(A, TrivialWritableVectorView<TX>(X));
}
template <class TA, class TX>
	typename IsWritableMatrix<typename TA::writable_matrix,
	typename IsWritableVector<typename TX::writable_vector,
MatrixOpStatus
	>::type>::type
SolveDestructive(TA &A, TX &X){
	return SolveDestructive(TrivialWritableMatrixView<TA>(A), TrivialWritableVectorView<TX>(X));
}


//// Solve(A, B, X) - A*X == B

//// Invert(A)

template <class TA, class TAinv>
	typename IsWritableMatrixView<typename TA::writable_matrix,
	typename IsWritableMatrixView<typename TAinv::writable_matrix,
MatrixOpStatus
	>::type>::type
InvertDestructive(const TA &A, const TAinv &Ainv){
	Fill(Ainv, typename TAinv::value_type(0));
	Fill(Diagonal(Ainv), typename TAinv::value_type(1));
	return SolveDestructive(A, Ainv);
}
template <class TA, class TAinv>
	typename IsWritableMatrix<typename TA::writable_matrix,
	typename IsWritableMatrixView<typename TAinv::writable_matrix,
MatrixOpStatus
	>::type>::type
InvertDestructive(TA &A, const TAinv &Ainv){
	return InvertDestructive(TrivialWritableMatrixView<TA>(A), Ainv);
}
template <class TA, class TAinv>
	typename IsWritableMatrixView<typename TA::writable_matrix,
	typename IsWritableMatrix<typename TAinv::writable_matrix,
MatrixOpStatus
	>::type>::type
InvertDestructive(const TA &A, TAinv &Ainv){
	return InvertDestructive(A, TrivialWritableMatrixView<TAinv>(Ainv));
}
template <class TA, class TAinv>
	typename IsWritableMatrix<typename TA::writable_matrix,
	typename IsWritableMatrix<typename TAinv::writable_matrix,
MatrixOpStatus
	>::type>::type
InvertDestructive(TA &A, TAinv &Ainv){
	return InvertDestructive(TrivialWritableMatrixView<TA>(A), TrivialWritableMatrixView<TAinv>(Ainv));
}

//// Eigensystem(A, L, X) - A*X == X*diag(L)

template <class TA, class TAlloc>
MatrixOpStatus Eigensystem(const TMatrix<TA,TAlloc> &A, TVector<TA,TAlloc> &Eval, TMatrix<TA,TAlloc> &Evec){
#ifdef USE_MATRIX_OPS_CHECKS
	if(A.Rows() != A.Cols() || A.Rows() != Eval.size() || Evec.Rows() != A.Rows() || Evec.Cols() != A.Cols()){ return DIMENSION_MISMATCH; }
#endif
#ifdef USE_MATRIX_OPS_ASSERTS
	assert(A.Rows() == A.Cols());
	assert(A.Rows() == Evec.Rows());
	assert(A.Cols() == Evec.Cols());
	assert(A.Rows() == Eval.size());
#endif
	return OK;
}

#endif // USE_ADVANCED_MATRIX_OPS



}; // namespace MatrixOps

#endif // _MATRIX_OPS_H_
