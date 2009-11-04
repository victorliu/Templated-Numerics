#ifndef _MATRIX_VIEW_H_
#define _MATRIX_VIEW_H_

#include <cstdlib>
#include "TVector.h"
#include "TMatrix.h"

// Preprocessor flags:
//   USE_MATRIX_ASSERTS
//   USE_COMPLEX_MATRICES

#ifdef USE_MATRIX_ASSERTS
# include <cassert>
#endif

#ifdef USE_COMPLEX_MATRICES
#include <complex>
#endif

#ifdef USE_COMPLEX_MATRICES

template <class ViewClass>
class ConjugateView : public VectorViewBase<typename ViewClass::value_type>{
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	ConjugateView(ViewClass V):view(V){}
	value_type  operator[](size_t row) const{ return std::conj(view[row]); }
	value_type& operator[](size_t row){ // illegal to use
		static value_type black_hole;
#ifdef USE_MATRIX_ASSERTS
		assert(0);
#endif
		return black_hole;
	}
	size_t size() const{ return view.size(); }
};

#endif // USE_COMPLEX_MATRICES

template <class ViewClass>
class DiagonalView : public VectorViewBase<typename ViewClass::value_type>{
protected:
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	
	DiagonalView(ViewClass V):view(V){}
	value_type& operator[](size_t row){ return view(row,row); }
	value_type  operator[](size_t row) const{ return view(row,row); }
	size_t size() const{ return ((view.Rows() < view.Cols()) ? view.Rows() : view.Cols()); }
};
template <class MatrixType>
DiagonalView<TrivialMatrixView<MatrixType> > Diagonal(MatrixType &M){
	return DiagonalView<TrivialMatrixView<MatrixType> >(TrivialMatrixView<MatrixType>(M));
}
template <class ViewType>
DiagonalView<ViewType> Transpose(const ViewType &V){
	return DiagonalView<ViewType>(V);
}

template <class ViewType>
class TransposeView : public MatrixViewBase<typename ViewType::value_type>{
	ViewType view;
public:
	typedef typename ViewType::value_type value_type;
	typedef ViewType parent_view;

	TransposeView(ViewType View):view(View){}
	value_type  operator()(size_t row, size_t col) const{ return view(col, row); }
	value_type& operator()(size_t row, size_t col)      { return view(col, row); }
	size_t Rows() const{ return view.Cols(); }
	size_t Cols() const{ return view.Rows(); }
};
template <class MatrixType>
TransposeView<TrivialMatrixView<MatrixType> > Transpose(MatrixType &M){
	return TransposeView<TrivialMatrixView<MatrixType> >(TrivialMatrixView<MatrixType>(M));
}
template <class ViewType>
TransposeView<ViewType> Transpose(const ViewType &V){
	return TransposeView<ViewType>(V);
}

#ifdef USE_COMPLEX_MATRICES

template <class ViewClass>
class ConjugateTransposeView : public MatrixViewBase<typename ViewClass::value_type>{
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	
	ConjugateTransposeView(ViewClass V):view(V){}
	value_type  operator()(size_t row, size_t col) const{ return std::conj(view(col, row)); }
	value_type& operator()(size_t row, size_t col){ // illegal to use
		static value_type black_hole;
#ifdef USE_MATRIX_ASSERTS
		assert(0);
#endif
		return black_hole;
	}
	size_t Rows() const{ return view.Cols(); }
	size_t Cols() const{ return view.Rows(); }
};
template <class MatrixType>
ConjugateTransposeView<TrivialMatrixView<MatrixType> > ConjugateTranspose(MatrixType &M){
	return ConjugateTransposeView<TrivialMatrixView<MatrixType> >(TrivialMatrixView<MatrixType>(M));
}
template <class ViewType>
ConjugateTransposeView<ViewType> ConjugateTranspose(const ViewType &V){
	return ConjugateTransposeView<ViewType>(V);
}

#endif // USE_COMPLEX_MATRICES

#endif // _MATRIX_VIEW_H_
