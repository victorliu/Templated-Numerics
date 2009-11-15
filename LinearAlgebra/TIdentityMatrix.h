#ifndef _TIDENTITY_MATRIX_H_
#define _TIDENTITY_MATRIX_H_

#include "TMatrix.h"
#include "TVector.h"

template <typename NumericType>
class TIdentityMatrix : public TMatrixBase<NumericType>{
	size_t rows;
	TIdentityMatrix(const TIdentityMatrix &M){}
	TIdentityMatrix& operator=(const TIdentityMatrix &M){ return *this; }
public:
	typedef NumericType value_type;
	
	TIdentityMatrix():rows(0){}
	TIdentityMatrix(size_t r):rows(r){}
	virtual ~TIdentityMatrix(){}
	
	void Resize(size_t nRows){ rows = nRows; }
	
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return rows; }
	value_type  operator[](size_t row) const{ return value_type(1); }
	value_type& operator[](size_t row){ return value_type(1); }
	value_type& operator()(size_t row, size_t col){
		static value_type one(1);
		if(row == col){ return one; }
		static value_type black_hole;
#ifdef USE_MATRIX_ASSERTS
		assert(0);
#endif
		return black_hole;
	}
	value_type operator()(size_t row, size_t col) const{ return (row == col) ? value_type(1) : value_type(0); }
};

#endif // _TIDENTITY_MATRIX_H_
