#ifndef _TDIAGONAL_MATRIX_H_
#define _TDIAGONAL_MATRIX_H_

#include "TMatrix.h"
#include "TVector.h"

template <typename NumericType, class TAllocator = std::allocator<NumericType> >
class TDiagonalMatrix : public TMatrixBase<NumericType>{
	NumericType *v;
	size_t rows;
	TAllocator allocator;
public:
	typedef NumericType value_type;
	
	TDiagonalMatrix():v(NULL),rows(0){}
	TDiagonalMatrix(size_t r):v(NULL),rows(r){
		v = allocator.allocate(r);
	}
	TDiagonalMatrix(size_t r, const value_type &init_val):v(NULL),rows(r){
		v = allocator.allocate(r);
		std::uninitialized_fill_n(v, r, init_val);
	}
	TDiagonalMatrix(const TDiagonalMatrix &M):v(NULL),rows(M.Rows()){
		v = allocator.allocate(rows);
		std::uninitialized_copy(M.Raw(), M.Raw()+rows, v);
	}
	TDiagonalMatrix& operator=(const TDiagonalMatrix &M){
		if(this != &M){
			Resize(M.Rows()); // cannot update rows and cols yet
			std::uninitialized_copy(M.Raw(), M.Raw()+rows, v);
		}
		return *this;
	}
	// Extract diagonal
	template <class T>
	TDiagonalMatrix(const TVectorBase<T> &V):v(NULL),rows(V.size()){
		v = allocator.allocate(rows);
		for(size_t i = 0; i < rows; ++i){
			(*this)[i] = V[i];
		}
	}
	virtual ~TDiagonalMatrix(){
		if(NULL != v){ allocator.deallocate(v, rows); }
	}
	
	void Resize(size_t nRows){
		if(NULL != v){ allocator.deallocate(v, rows); }
		rows = nRows;
		v = allocator.allocate(rows);
	}
	
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return rows; }
	value_type  operator[](size_t row) const{ return v[row]; }
	value_type& operator[](size_t row){ return v[row]; }
	value_type& operator()(size_t row, size_t col){
		static value_type black_hole;
		if(row == col){ return v[row]; }
#ifdef USE_MATRIX_ASSERTS
		assert(0);
#endif
		return black_hole;
	}
	value_type operator()(size_t row, size_t col) const{ return (row == col) ? v[row] : value_type(0); }
	
	void Invert(){
		for(size_t i = 0; i < rows; ++i){
			v[i] = value_type(1)/v[i];
		}
	}
	
	// Raw interface
	value_type* Raw() const{ return v; }
};

#endif // _TDIAGONAL_MATRIX_H_
