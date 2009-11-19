#ifndef _MATRIX_VIEWS_H_
#define _MATRIX_VIEWS_H_

#include <cstddef>
#include "MatrixInterfaces.h"

// Preprocessor flags:
//   USE_COMPLEX_MATRICES

#include <cassert>

#ifdef USE_COMPLEX_MATRICES
#include <complex>
#endif

// The trivial views provide adapters from matrices to matrix views,
// and correspondingly for vectors.
template <class MatrixType>
class TrivialReadableMatrixView : public ReadableMatrix<typename MatrixType::value_type>
{
	const MatrixType &M;
public:
	typedef typename MatrixType::value_type value_type;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	
	TrivialReadableMatrixView(const MatrixType &mat):M(mat){}
	value_type operator()(size_t row, size_t col) const{ return M(row, col); }
	size_t Rows() const{ return M.Rows(); }
	size_t Cols() const{ return M.Cols(); }
};
template <class MatrixType>
class TrivialWritableMatrixView : public WritableMatrixView<typename MatrixType::value_type>
{
	MatrixType &M;
public:
	typedef TrivialWritableMatrixView<MatrixType> self;
	typedef typename MatrixType::value_type value_type;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	typedef WritableMatrixView<value_type> writable_matrix;
	
	TrivialWritableMatrixView(MatrixType &mat):M(mat){}
	value_type& operator()(size_t row, size_t col) const{ return M(row,col); }
	size_t Rows() const{ return M.Rows(); }
	size_t Cols() const{ return M.Cols(); }
};

template <class VectorType>
class TrivialReadableVectorView : public ReadableVector<typename VectorType::value_type>
{
	const VectorType &V;
public:
	typedef typename VectorType::value_type value_type;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	
	TrivialReadableVectorView(const typename VectorType::value_type &vec):V(vec){}
	value_type operator[](size_t row) const{ return V[row]; }
	size_t size() const{ return V.size(); }
};
template <class VectorType>
class TrivialWritableVectorView : public WritableVectorView<typename VectorType::value_type>
{
	VectorType &V;
public:
	typedef typename VectorType::value_type value_type;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVectorView<value_type> writable_vector;
	
	TrivialWritableVectorView(VectorType &vec):V(vec){}
	value_type& operator[](size_t row) const{ return V[row]; }
	size_t size() const{ return V.size(); }
};

// Represents some constant scaling of all elements
template <class ViewClass>
class ScaledMatrixView : public ReadableMatrix<typename ViewClass::value_type>
{
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
protected:
	ViewClass view;
	value_type c;
public:
	ScaledMatrixView(parent_view V, const value_type &scale):view(V), c(scale){}
	value_type operator()(size_t row, size_t col) const{ return c*view(row, col); }
	size_t Rows() const{ return view.Rows(); }
	size_t Cols() const{ return view.Cols(); }
};
template <class ViewType>
	typename IsMatrixView<typename ViewType::view_type,
ScaledMatrixView<ViewType>
	>::type
Scaled(const ViewType &M, const typename ViewType::value_type &scale){
	return ScaledMatrixView<ViewType>(M, scale);
}
template <class MatrixType>
	typename IsMatrixNonView<typename MatrixType::non_view_type,
ScaledMatrixView<TrivialReadableMatrixView<MatrixType> >
	>::type
Scaled(const MatrixType &M, const typename MatrixType::value_type &scale){
	return ScaledMatrixView<TrivialReadableMatrixView<MatrixType> >(TrivialReadableMatrixView<MatrixType>(M), scale);
}

template <class ViewClass>
class ScaledVectorView : public ReadableVector<typename ViewClass::value_type>
{
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
protected:
	ViewClass view;
	value_type c;
public:
	ScaledVectorView(parent_view V, const value_type &scale):view(V), c(scale){}
	value_type operator[](size_t row) const{ return c*view[row]; }
	size_t size() const{ return view.size(); }
};
template <class ViewType>
	typename IsVectorView<typename ViewType::view_type,
ScaledVectorView<ViewType>
	>::type
Scaled(const ViewType &M, const typename ViewType::value_type &scale){
	return ScaledVectorView<ViewType>(M, scale);
}
template <class VectorType>
	typename IsVectorNonView<typename VectorType::non_view_type,
ScaledVectorView<TrivialReadableVectorView<VectorType> >
	>::type
Scaled(const VectorType &M, const typename VectorType::value_type &scale){
	return ScaledVectorView<TrivialReadableVectorView<VectorType> >(TrivialReadableVectorView<VectorType>(M), scale);
}


// The Sub views select portions of matrices/vectors
template <class ViewClass>
class SubMatrixView : public WritableMatrixView<typename ViewClass::value_type>
{
	ViewClass view;
	size_t row_start, col_start, rows, cols;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	typedef WritableMatrixView<value_type> writable_matrix;
	
	SubMatrixView(parent_view V, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols):view(V),row_start(RowStart),col_start(ColStart),rows(nRows),cols(nCols){}
	value_type& operator()(size_t row, size_t col) const{ return view(row+row_start, col+col_start); }
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
};
template <class ViewType>
	typename IsMatrixView<typename ViewType::view_type,
SubMatrixView<ViewType>
	>::type
SubMatrix(const ViewType &M, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols){
	return SubMatrixView<ViewType>(M, RowStart, ColStart, nRows, nCols);
}
template <class MatrixType>
	typename IsMatrixNonView<typename MatrixType::non_view_type,
SubMatrixView<TrivialWritableMatrixView<MatrixType> >
	>::type
SubMatrix(MatrixType &M, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols){
	return SubMatrixView<TrivialWritableMatrixView<MatrixType> >(TrivialWritableMatrixView<MatrixType>(M), RowStart, ColStart, nRows, nCols);
}

template <class ViewClass>
class SubVectorView : public WritableVectorView<typename ViewClass::value_type>
{
	ViewClass view;
	size_t row_start, rows;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVectorView<value_type> writable_vector;
	
	SubVectorView(parent_view V, size_t RowStart, size_t nRows):view(V),row_start(RowStart),rows(nRows){}
	value_type& operator[](size_t row) const{ return view[row+row_start]; }
	size_t size() const{ return rows; }
};
template <class ViewType>
	typename IsVectorView<typename ViewType::view_type,
SubVectorView<ViewType>
	>::type
SubVector(const ViewType &V, size_t RowStart, size_t Len){
	return SubVectorView<ViewType>(V, RowStart, Len);
}
template <class VectorType>
	typename IsVectorNonView<typename VectorType::non_view_type,
SubVectorView<TrivialWritableVectorView<VectorType> >
	>::type
SubVector(VectorType &V, size_t RowStart, size_t Len){
	return SubVectorView<TrivialWritableVectorView<VectorType> >(TrivialWritableVectorView<VectorType>(V), RowStart, Len);
}

template <class ViewClass>
class DiagonalView : public WritableVectorView<typename ViewClass::value_type>{
protected:
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVectorView<value_type> writable_vector;
	
	DiagonalView(ViewClass V):view(V){}
	value_type& operator[](size_t row) const{ return view(row,row); }
	size_t size() const{ return ((view.Rows() < view.Cols()) ? view.Rows() : view.Cols()); }
};
template <class ViewType>
	typename IsMatrixView<typename ViewType::view_type,
DiagonalView<ViewType>
	>::type
Diagonal(const ViewType &V){
	return DiagonalView<ViewType>(V);
}
template <class MatrixType>
	typename IsMatrixNonView<typename MatrixType::non_view_type,
	typename IsWritableMatrix<typename MatrixType::writable_matrix,
DiagonalView<TrivialWritableMatrixView<MatrixType> >
	>::type>::type
Diagonal(MatrixType &V){
	return DiagonalView<TrivialWritableMatrixView<MatrixType> >(TrivialWritableMatrixView<MatrixType>(V));
}


template <class ViewClass>
class ColumnView : public WritableVectorView<typename ViewClass::value_type>{
	ViewClass view;
	size_t col;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVectorView<value_type> writable_vector;
	
	ColumnView(ViewClass V, size_t nCol):view(V),col(nCol){}
	value_type& operator[](size_t row) const{ return view(row, col); }
	size_t size() const{ return view.Rows(); }
};
template <class ViewType>
	typename IsMatrixView<typename ViewType::view_type,
ColumnView<ViewType>
	>::type
GetColumn(const ViewType &V, size_t col){
	return ColumnView<ViewType>(V, col);
}
template <class MatrixType>
	typename IsMatrixNonView<typename MatrixType::non_view_type,
ColumnView<TrivialWritableMatrixView<MatrixType> >
	>::type
GetColumn(MatrixType &M, size_t col){
	return ColumnView<TrivialWritableMatrixView<MatrixType> >(TrivialWritableMatrixView<MatrixType>(M), col);
}

template <class ViewClass>
class RowView : public WritableVectorView<typename ViewClass::value_type>{
	ViewClass view;
	size_t row;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVectorView<value_type> writable_vector;
	
	RowView(ViewClass V, size_t nRow):view(V),row(nRow){}
	value_type& operator[](size_t col) const{ return view(row, col); }
	size_t size() const{ return view.Cols(); }
};
template <class ViewType>
	typename IsMatrixView<typename ViewType::view_type,
RowView<ViewType>
	>::type
GetRow(const ViewType &V, size_t row){
	return RowView<ViewType>(V, row);
}
template <class MatrixType>
	typename IsMatrixNonView<typename MatrixType::non_view_type,
RowView<TrivialWritableMatrixView<MatrixType> >
	>::type
GetRow(MatrixType &M, size_t row){
	return RowView<TrivialWritableMatrixView<MatrixType> >(TrivialWritableMatrixView<MatrixType>(M), row);
}



template <class ViewClass>
class TransposeView : public WritableMatrixView<typename ViewClass::value_type>{
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	typedef WritableMatrixView<value_type> writable_matrix;

	TransposeView(ViewClass View):view(View){}
	value_type& operator()(size_t row, size_t col) const{ return view(col, row); }
	size_t Rows() const{ return view.Cols(); }
	size_t Cols() const{ return view.Rows(); }
};
template <class ViewType>
	typename IsMatrixView<typename ViewType::view_type,
TransposeView<ViewType>
	>::type
Transpose(const ViewType &V){
	return TransposeView<ViewType>(V);
}
template <class MatrixType>
	typename IsMatrixNonView<typename MatrixType::non_view_type,
TransposeView<TrivialWritableMatrixView<MatrixType> >
	>::type
Transpose(MatrixType &M){
	return TransposeView<TrivialWritableMatrixView<MatrixType> >(TrivialWritableMatrixView<MatrixType>(M));
}


#ifdef USE_COMPLEX_MATRICES

template <class ViewClass>
class ConjugateVectorView : public ReadableVector<typename ViewClass::value_type>{
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;

	ConjugateVectorView(ViewClass V):view(V){}
	value_type  operator[](size_t row) const{ return std::conj(view[row]); }
	size_t size() const{ return view.size(); }
};
template <class ViewType>
	typename IsMatrixView<typename ViewType::view_type,
ConjugateVectorView<ViewType>
	>::type
Conjugate(const ViewType &V){
	return ConjugateVectorView<ViewType>(V);
}
template <class VectorType>
	typename IsMatrixNonView<typename VectorType::non_view_type,
ConjugateVectorView<TrivialWritableMatrixView<VectorType> >
	>::type
Conjugate(const VectorType &V){
	return ConjugateVectorView<TrivialReadableVectorView<VectorType> >(TrivialReadableVectorView<VectorType>(V));
}


template <class ViewClass>
class ConjugateTransposeView : public ReadableMatrix<typename ViewClass::value_type>{
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	
	ConjugateTransposeView(ViewClass V):view(V){}
	value_type  operator()(size_t row, size_t col) const{ return std::conj(view(col, row)); }
	size_t Rows() const{ return view.Cols(); }
	size_t Cols() const{ return view.Rows(); }
};
template <class ViewType>
	typename IsMatrixView<typename ViewType::view_type,
ConjugateTransposeView<ViewType>
	>::type
ConjugateTranspose(const ViewType &V){
	return ConjugateTransposeView<ViewType>(V);
}
template <class MatrixType>
	typename IsMatrixNonView<typename MatrixType::non_view_type,
ConjugateTransposeView<TrivialWritableMatrixView<MatrixType> >
	>::type
ConjugateTranspose(const MatrixType &M){
	return ConjugateTransposeView<TrivialReadableMatrixView<MatrixType> >(TrivialReadableMatrixView<MatrixType>(M));
}

#endif // USE_COMPLEX_MATRICES



#endif // _MATRIX_VIEWS_H_
