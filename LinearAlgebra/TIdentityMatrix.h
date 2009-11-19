#ifndef _TIDENTITY_MATRIX_H_
#define _TIDENTITY_MATRIX_H_

#define USING_TIDENTITY_MATRIX

#include "MatrixInterfaces.h"

template <typename NumericType>
class TIdentityMatrix : public ReadableMatrix<NumericType>{
	size_t rows;
	TIdentityMatrix& operator=(const TIdentityMatrix &M){ return *this; }
public:
	typedef NumericType value_type;
	typedef ReadableMatrix<value_type> non_view_type;
	typedef ReadableMatrix<value_type> readable_matrix;
	typedef WritableMatrix<value_type> non_writable_matrix;
	
	TIdentityMatrix():rows(0){}
	TIdentityMatrix(size_t r):rows(r){}
	TIdentityMatrix(const TIdentityMatrix &M):rows(M.rows){}
	virtual ~TIdentityMatrix(){}
	
	void Resize(size_t nRows){ rows = nRows; }
	
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return rows; }
	value_type  operator[](size_t row) const{ return value_type(1); }
	value_type operator()(size_t row, size_t col) const{ return (row == col) ? value_type(1) : value_type(0); }
};

#endif // _TIDENTITY_MATRIX_H_
