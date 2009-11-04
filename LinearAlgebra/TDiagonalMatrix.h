#ifndef _TDIAGONAL_MATRIX_H_
#define _TDIAGONAL_MATRIX_H_

#include "MatrixView.h"


// View for dense matrices
template <class T>
class TDiagonalMatrixView : public MatrixViewBase<T>{
protected:
	T* A;
	size_t rows;
public:
	typedef T value_type;
	typedef TDiagonalMatrixView<value_type> matrix_type;
	
	TDiagonalMatrixView(value_type* DataPtr, size_t nRows):
		A(DataPtr),
		rows(nRows){
	}
	value_type& operator()(size_t row, size_t col)      { return A[row]; }
	value_type  operator()(size_t row, size_t col) const{ return (row == col) ? A[row] : value_type(0); }
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return rows; }
	
	size_t LeadingDimension() const{ return rows; }
};

template <typename NumericType, class TAllocator = std::allocator<NumericType> >
class TDiagonalMatrix : public TMatrixBase<NumericType>{
	NumericType *v;
	size_t rows;
	TAllocator allocator;
public:
	typedef NumericType value_type;
	typedef TDiagonalMatrixView<value_type> View;
	typedef TDiagonalMatrix<value_type> matrix_type;
	
	TDiagonalMatrix():v(NULL),rows(0){}
	TDiagonalMatrix(size_t r):v(NULL),rows(r){
		v = allocator.allocate(r);
	}
	TDiagonalMatrix(size_t r, const value_type &init_val):v(NULL),rows(r){
		v = allocator.allocate(r);
		std::uninitialized_fill_n(v, r, init_val);
	}
	TDiagonalMatrix(const TDiagonalMatrix &V):v(NULL),rows(M.Rows()){
		v = allocator.allocate(rows);
		std::uninitialized_copy(M.Raw(), M.Raw()+rows, v);
	}
	TDiagonalMatrix& operator=(const TDiagonalMatrix &M){
		if(this != &M){
			Resize(M.size()); // cannot update rows and cols yet
			std::uninitialized_copy(M.Raw(), M.Raw()+rows, v);
		}
		return *this;
	}
	TDiagonalMatrix(const TDiagonalMatrixView<value_type> &M):v(NULL),rows(M.Rows()){
		v = allocator.allocate(rows);
		std::uninitialized_copy(M.Raw(), M.Raw()+rows, v);
	}
	// Extract diagonal
	template <class MatrixLikeType>
	TDiagonalMatrix(const MatrixLikeType &M):v(NULL),rows(V.Rows()){
		if(rows > M.Cols()){ rows = M.Cols(); }
		v = allocator.allocate(rows);
		for(size_t i = 0; i < rows; ++i){
			(*this)[i] = M(i,i);
		}
	}
	~TDiagonalMatrix(){
		allocator.deallocate(v, rows);
	}
	
	void Resize(size_t nRows){
		allocator.deallocate(v, rows);
		rows = nRows;
		v = allocator.allocate(rows);
	}
	
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return rows; }
	const NumericType& operator[](size_t row) const{ return v[row]; }
	NumericType& operator[](size_t row){ return v[row]; }
	value_type& operator()(size_t row, size_t col)      { return v[row]; }
	const value_type& operator()(size_t row, size_t col) const{
		static const value_type zero(0);
		return (row == col) ? v[row] : zero;
	}
	operator TDiagonalMatrixView<value_type>(){ return TDiagonalMatrixView<value_type>(v, rows); }
	
	// Raw interface
	NumericType* Raw() const{ return v; }
};

#endif // _TDIAGONAL_MATRIX_H_
