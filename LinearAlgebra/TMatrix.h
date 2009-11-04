#ifndef _TMATRIX_H_
#define _TMATRIX_H_

#include "MatrixView.h"

template <typename NumericType>
class TMatrixBase{
public:
	typedef NumericType value_type;
	typedef TMatrixBase<value_type> matrix_type;
	// Implementers must define this:
	//typedef MatrixViewBase<T> View;
	
	virtual size_t Rows() const = 0;
	virtual size_t Cols() const = 0;
	virtual const NumericType& operator()(size_t row, size_t col) const = 0;
	virtual       NumericType& operator()(size_t row, size_t col)       = 0;
	
	//virtual operator MatrixView<value_type>&() = 0;
};



// View for dense matrices
template <class T>
class TMatrixView : public MatrixViewBase<T>{
protected:
	T* A;
	size_t rows, cols, col_stride;
public:
	typedef T value_type;
	typedef TMatrixView<value_type> matrix_type;
	
	TMatrixView(value_type* DataPtr, size_t nRows, size_t nCols, size_t lda):
		A(DataPtr),
		rows(nRows),
		cols(nCols),
		col_stride(lda){
	}
	value_type& operator()(size_t row, size_t col)      { return A[col*col_stride+row]; }
	value_type  operator()(size_t row, size_t col) const{ return A[col*col_stride+row]; }
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
	
	size_t LeadingDimension() const{ return col_stride; }
	
	friend class SubMatrixView<TMatrixView>;
};

// Specialization of submatrix view for dense matrices
template <typename T>
class SubMatrixView<TMatrixView<T> > : public MatrixViewBase<typename TMatrixView<T>::value_type>{
	TMatrixView view;
public:
	typedef typename TMatrixView<T>::value_type value_type;
	typedef SubMatrixView<value_type> matrix_type;
	
	SubMatrixView(const TMatrixView<T> &view, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols):TMatrixView<T>(&view(RowStart,ColStart), nRows, nCols, view.col_stride){}
	value_type  operator()(size_t row, size_t col) const{ return view(row, col); }
	value_type& operator()(size_t row, size_t col)      { return view(row, col); }
	size_t Rows() const{ return view.Rows(); }
	size_t Cols() const{ return view.Cols(); }
	operator TMatrixView<T>(){ return view; }
};


template <typename NumericType, class TAllocator = std::allocator<NumericType> >
class TMatrix : public TMatrixBase<NumericType>{
	NumericType *A;
	size_t rows, cols;
	TAllocator allocator;
public:
	typedef NumericType value_type;
	typedef TMatrixView<value_type> View;
	typedef TMatrix<value_type> matrix_type;
	
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
		if(this != &M){
			Resize(M.Rows(), M.Cols()); // cannot update rows and cols yet
			std::uninitialized_copy(M.Raw(), M.Raw()+rows*cols, A);
		}
		return *this;
	}
	TMatrix(const TMatrixView<value_type> &M):A(NULL),rows(M.Rows()),cols(M.Cols()){
		A = allocator.allocate(rows*cols);
		std::uninitialized_copy(M.Raw(), M.Raw()+rows*cols, A);
	}
	template <class MatrixLikeType>
	TMatrix(const MatrixLikeType &M):A(NULL),rows(M.Rows()),cols(M.Cols()){
		A = allocator.allocate(rows*cols);
		for(size_t j = 0; j < cols; ++j){
			for(size_t i = 0; i < rows; ++i){
				(*this)(i,j) = M(i,j);
			}
		}
	}
	~TMatrix(){
		allocator.deallocate(A, rows*cols);
	}
	
	void Resize(size_t nRows, size_t nCols){
		allocator.deallocate(A, rows*cols);
		rows = nRows; cols = nCols;
		A = allocator.allocate(rows*cols);
	}
	
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
	const NumericType& operator()(size_t row, size_t col) const{
		return A[rows*col + row];
	}
	NumericType& operator()(size_t row, size_t col){
		return A[rows*col + row];
	}
	operator TMatrixView<value_type>(){
		return TMatrixView<value_type>(A, rows, cols, rows);
	}
	TMatrixView<value_type> SubMatrix(size_t starting_row, size_t starting_col, size_t num_rows, size_t num_cols){
		return TMatrixView<value_type>(&(A(starting_row, starting_col)), num_rows, num_cols, rows);
	}
	
	// Raw interface
	NumericType* Raw() const{ return A; }
	size_t LeadingDimension() const{ return rows; }
};

#endif // _TMATRIX_H_
