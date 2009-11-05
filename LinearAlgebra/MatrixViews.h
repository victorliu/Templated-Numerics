#ifndef _MATRIX_VIEWS_H_
#define _MATRIX_VIEWS_H_

#define USING_MATRIX_VIEWS

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



template <class ViewClass>
class SubVectorView : public VectorViewBase<typename ViewClass::value_type>{
	ViewClass view;
	size_t row_start, rows;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	
	SubVectorView(ViewClass V, size_t RowStart, size_t nRows):view(V),row_start(RowStart),rows(nRows){}
	value_type  operator[](size_t row) const{ return view[row+row_start]; }
	value_type& operator[](size_t row)      { return view[row+row_start]; }
	size_t size() const{ return rows; }
};
template <class VectorType>
SubVectorView<TrivialVectorView<VectorType> > SubVector(VectorType &V, size_t start){
	return SubVectorView<TrivialVectorView<VectorType> >(TrivialVectorView<VectorType>(V), start, V.size()-start);
}
template <class ViewType>
SubVectorView<ViewType> SubVector(const ViewType &V, size_t start){
	return SubVectorView<ViewType>(V, start, V.size()-start);
}
template <class VectorType>
SubVectorView<TrivialVectorView<VectorType> > SubVector(VectorType &V, size_t start, size_t len){
	return SubVectorView<TrivialVectorView<VectorType> >(TrivialVectorView<VectorType>(V), start, len);
}
template <class ViewType>
SubVectorView<ViewType> SubVector(const ViewType &V, size_t start, size_t len){
	return SubVectorView<ViewType>(V, start, len);
}

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



template <class ViewClass>
class SubMatrixView : public MatrixViewBase<typename ViewClass::value_type>{
	ViewClass view;
	size_t row_start, col_start, rows, cols;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	
	SubMatrixView(ViewClass V, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols):view(V),row_start(RowStart),col_start(ColStart),rows(nRows),cols(nCols){}
	value_type  operator()(size_t row, size_t col) const{ return view(row+row_start, col+col_start); }
	value_type& operator()(size_t row, size_t col)      { return view(row+row_start, col+col_start); }
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
};
template <class MatrixType>
SubMatrixView<TrivialMatrixView<MatrixType> > SubMatrix(MatrixType &M, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols){
	return SubMatrixView<TrivialMatrixView<MatrixType> >(TrivialMatrixView<MatrixType>(M), RowStart, ColStart, nRows, nCols);
}
template <class ViewType>
SubMatrixView<ViewType> SubMatrix(const ViewType &M, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols){
	return SubMatrixView<ViewType>(M, RowStart, ColStart, nRows, nCols);
}

template <class ViewClass>
class ColumnView : public VectorViewBase<typename ViewClass::value_type>{
	ViewClass view;
	size_t col;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	
	ColumnView(ViewClass V, size_t nCol):view(V),col(nCol){}
	value_type  operator[](size_t row) const{ return view(row, col); }
	value_type& operator[](size_t row)      { return view(row, col); }
	size_t size() const{ return view.Rows(); }
};
template <class MatrixType>
ColumnView<TrivialMatrixView<MatrixType> > GetColumn(MatrixType &M, size_t col){
	return ColumnView<TrivialMatrixView<MatrixType> >(TrivialMatrixView<MatrixType>(M), col);
}
template <class ViewType>
ColumnView<ViewType> GetColumn(const ViewType &V, size_t col){
	return ColumnView<ViewType>(V, col);
}

template <class ViewClass>
class RowView : public VectorViewBase<typename ViewClass::value_type>{
	ViewClass view;
	size_t row;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	
	RowView(ViewClass V, size_t nRow):view(V),row(nRow){}
	value_type  operator[](size_t col) const{ return view(row, col); }
	value_type& operator[](size_t col)      { return view(row, col); }
	size_t size() const{ return view.Cols(); }
};
template <class MatrixType>
RowView<TrivialMatrixView<MatrixType> > GetRow(MatrixType &M, size_t row){
	return RowView<TrivialMatrixView<MatrixType> >(TrivialMatrixView<MatrixType>(M), row);
}
template <class ViewType>
RowView<ViewType> GetRow(const ViewType &V, size_t row){
	return RowView<ViewType>(V, row);
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

#endif // _MATRIX_VIEWS_H_
