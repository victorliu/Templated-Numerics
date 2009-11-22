#ifndef _TMATRIX_H_
#define _TMATRIX_H_

// Preprocessor flags:

#include "MatrixInterfaces.h"
#include <memory>
#include <cassert>
#include "MatrixViews.h"

template <typename NumericType, class TAllocator = std::allocator<NumericType> >
class TMatrix : public WritableMatrix<NumericType>{
public:
	typedef NumericType value_type;
	typedef TMatrix<value_type,TAllocator> self;
	
	typedef ReadableMatrix<value_type> non_view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	typedef WritableMatrix<value_type> writable_matrix;
protected:
	NumericType *A;
	size_t rows, cols;
	TAllocator allocator;
public:
	TMatrix():A(NULL),rows(0),cols(0){}
	TMatrix(size_t r, size_t c):A(NULL),rows(r),cols(c){
		A = allocator.allocate(r*c);
	}
	TMatrix(size_t r, size_t c, const value_type &init_val):A(NULL),rows(r),cols(c){
		A = allocator.allocate(r*c);
		std::uninitialized_fill_n(A, r*c, init_val);
	}
	TMatrix(const TMatrix &M):A(NULL),rows(M.Rows()),cols(M.Cols()){
		A = allocator.allocate(rows*cols);
		std::uninitialized_copy(M.Raw(), M.Raw()+rows*cols, A);
	}
	TMatrix& operator=(const TMatrix &M){
		Resize(M.Rows(), M.Cols());
		for(size_t j = 0; j < cols; ++j){
			for(size_t i = 0; i < rows; ++i){
				(*this)(i,j) = M(i,j);
			}
		}
		return *this;
	}
	template <class MatrixLikeType>
	TMatrix& operator=(const MatrixLikeType &M){
		Resize(M.Rows(), M.Cols());
		for(size_t j = 0; j < cols; ++j){
			for(size_t i = 0; i < rows; ++i){
				(*this)(i,j) = M(i,j);
			}
		}
		return *this;
	}
	virtual ~TMatrix(){
		if(NULL != A){ allocator.deallocate(A, rows*cols); }
	}
	
	void Resize(size_t nRows, size_t nCols){
		if(NULL != A){ allocator.deallocate(A, rows*cols); }
		rows = nRows; cols = nCols;
		A = allocator.allocate(rows*cols);
	}
	
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
	NumericType operator()(size_t row, size_t col) const{
		assert(row < Rows());
		assert(col < Cols());
		return A[rows*col + row];
	}
	NumericType& operator()(size_t row, size_t col){
		assert(row < Rows());
		assert(col < Cols());
		return A[rows*col + row];
	}
	/*
	operator TMatrixView<value_type>(){
		return TMatrixView<value_type>(A, rows, cols, rows);
	}
	TMatrixView<value_type> SubMatrix(size_t starting_row, size_t starting_col, size_t num_rows, size_t num_cols){
		return TMatrixView<value_type>(&(A(starting_row, starting_col)), num_rows, num_cols, rows);
	}
	operator WritableMatrixArgument<self>(){
		return WritableMatrixArgument<self>(*this);
	}*/
	
	// Raw interface
	NumericType* Raw() const{ return A; }
	size_t LeadingDimension() const{ return rows; }
};




template <class T, class TAlloc>
class TrivialReadableMatrixView<TMatrix<T,TAlloc> > : public ReadableMatrix<T>
{
	const TMatrix<T,TAlloc> &M;
public:
	typedef T value_type;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	
	TrivialReadableMatrixView(const TMatrix<T,TAlloc> &mat):M(mat){}
	value_type operator()(size_t row, size_t col) const{ return M(row, col); }
	size_t Rows() const{ return M.Rows(); }
	size_t Cols() const{ return M.Cols(); }
	
	value_type *Raw() const{ return M.Raw(); }
	size_t LeadingDimension() const{ return M.LeadingDimension(); }
};
template <class T, class TAlloc>
class TrivialWritableMatrixView<TMatrix<T,TAlloc> > : public WritableMatrixView<T>
{
	TMatrix<T,TAlloc> &M;
public:
	typedef T value_type;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	typedef WritableMatrixView<value_type> writable_matrix;
	
	TrivialWritableMatrixView(TMatrix<T,TAlloc> &mat):M(mat){}
	void Set(size_t row, size_t col, const value_type &value) const{ M(row,col) = value; }
	const value_type& Get(size_t row, size_t col) const{ return M(row,col); }
	value_type& GetMutable(size_t row, size_t col) const{ return M(row,col); }
	size_t Rows() const{ return M.Rows(); }
	size_t Cols() const{ return M.Cols(); }
	
	value_type *Raw() const{ return M.Raw(); }
	size_t LeadingDimension() const{ return M.LeadingDimension(); }
};


// Specialization of submatrix view for dense matrices
template <typename T, class TAlloc>
class SubMatrixView<TrivialWritableMatrixView<TMatrix<T,TAlloc> > > : public WritableMatrixView<T>{
	T* A;
	size_t rows, cols, lda;
public:
	typedef T value_type;
	typedef TrivialWritableMatrixView<TMatrix<T,TAlloc> > parent_view;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	typedef WritableMatrixView<value_type> writable_matrix;
	
	SubMatrixView<TrivialWritableMatrixView<TMatrix<T,TAlloc> > >(TrivialWritableMatrixView<TMatrix<T,TAlloc> > View, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols):A(&View.GetMutable(RowStart,ColStart)),rows(nRows),cols(nCols),lda(View.LeadingDimension()){}
	void Set(size_t row, size_t col, const value_type &value) const{ A[row+col*lda] = value; }
	const value_type& Get(size_t row, size_t col) const{ return A[row+col*lda]; }
	value_type& GetMutable(size_t row, size_t col) const{ return A[row+col*lda]; }
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
	
	size_t LeadingDimension() const{ return lda; }
	value_type *Raw() const{ return A; }
};
// Specialization of submatrix view for dense matrices
template <typename T, class TAlloc>
class SubMatrixView<TrivialReadableMatrixView<TMatrix<T,TAlloc> > > : public ReadableMatrix<T>{
	T* A;
	size_t rows, cols, lda;
public:
	typedef T value_type;
	typedef TrivialWritableMatrixView<TMatrix<T,TAlloc> > parent_view;
	typedef ReadableMatrix<value_type> view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	typedef WritableMatrixView<value_type> writable_matrix;
	
	SubMatrixView<TrivialReadableMatrixView<TMatrix<T,TAlloc> > >(TrivialReadableMatrixView<TMatrix<T,TAlloc> > View, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols):A(&(View(RowStart,ColStart))),rows(nRows),cols(nCols),lda(View.LeadingDimension()){}
	value_type operator()(size_t row, size_t col) const{ return A[row+col*lda]; }
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
	
	size_t LeadingDimension() const{ return lda; }
	value_type *Raw() const{ return A; }
};
template <class T, class TAlloc>
SubMatrixView<TrivialWritableMatrixView<TMatrix<T,TAlloc> > > SubMatrix(TMatrix<T,TAlloc> &M, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols){
	return SubMatrixView<TrivialWritableMatrixView<TMatrix<T,TAlloc> > >(
		TrivialWritableMatrixView<TMatrix<T,TAlloc> >(M),
		RowStart, ColStart, nRows, nCols);
}
template <class T, class TAlloc>
SubMatrixView<TrivialReadableMatrixView<TMatrix<T,TAlloc> > > SubMatrix(const TMatrix<T,TAlloc> &M, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols){
	return SubMatrixView<TrivialReadableMatrixView<TMatrix<T,TAlloc> > >(
		TrivialReadableMatrixView<TMatrix<T,TAlloc> >(M),
		RowStart, ColStart, nRows, nCols);
}

#endif // _TMATRIX_H_
