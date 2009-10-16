#ifndef _TMATRIX_H_
#define _TMATRIX_H_

#include "MatrixView.h"

class MatrixBase{
};

template <typename NumericType>
class TMatrix{
	NumericType *A;
	size_t rows, cols;
	size_t lda;
public:
	Matrix();
	Matrix(size_t r, size_t c);
	Matrix(const Matrix &M);
	
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
	const NumericType& operator()(size_t row, size_t col) const{
		return A[lda*col + row];
	}
	NumericType& operator()(size_t row, size_t col){
		return A[lda*col + row];
	}
	operator MatrixViewBase(){
		return MatrixViewBase(A, rows, cols, lda);
	}
	MatrixViewBase SubMatrix(size_t starting_row, size_t starting_col, size_t num_rows, size_t num_cols){
		return MatrixViewBase(&(A(starting_row, starting_col), num_rows, num_cols, lda);
	}
	
	// Raw interface
	NumericType* Raw() const{ return A; }
};

#endif // _TMATRIX_H_
