#ifndef _TVECTOR_H_
#define _TVECTOR_H_

template <typename NumericType>
class TVectorBase{
public:
	typedef NumericType value_type;
	
	virtual size_t size() const = 0;
	size_t Rows() const{ return size(); }
	size_t Cols() const{ return 1; }
	value_type  operator()(size_t row, size_t col) const{
#ifdef USE_MATRIX_ASSERTS
		assert(1 == col);
#endif
		return (*this)[row];
	}
	value_type& operator()(size_t row, size_t col){
#ifdef USE_MATRIX_ASSERTS
		assert(1 == col);
#endif
		return (*this)[row];
	}
	virtual NumericType  operator[](size_t row) const = 0;
	virtual NumericType& operator[](size_t row)       = 0;
};


template <class T>
class VectorViewBase : public TVectorBase<T>{
public:
	typedef T value_type;
	//typedef ... parent_view;
};

template <class VectorType>
class TrivialVectorView : public VectorViewBase<typename VectorType::value_type>{
	VectorType &V;
public:
	typedef typename VectorType::value_type value_type;
	
	TrivialVectorView(VectorType &vec):V(vec){}
	value_type  operator[](size_t row) const{ return V[row]; }
	value_type& operator[](size_t row)      { return V[row]; }
	size_t size() const{ return V.size(); }
};

#include "MatrixViews.h"

// View for dense vectors
template <class T>
class TVectorView : public VectorViewBase<T>{
protected:
	T* V;
	size_t rows, stride;
public:
	typedef T value_type;
	
	TVectorView(T* DataPtr, size_t nRows, size_t Stride):
		V(DataPtr),
		rows(nRows),
		stride(Stride){
	}
	virtual T& operator[](size_t row)      { return V[stride*row]; }
	virtual T  operator[](size_t row) const{ return V[stride*row]; }
	virtual size_t size() const{ return rows; }
};

#include <memory>

template <typename NumericType, class TAllocator = std::allocator<NumericType> >
class TVector : public TVectorBase<NumericType>{
	NumericType *v;
	size_t rows;
	TAllocator allocator;
public:
	typedef NumericType value_type;
	
	TVector():v(NULL),rows(0){}
	TVector(size_t r):v(NULL),rows(r){
		v = allocator.allocate(r);
	}
	TVector(size_t r, const value_type &init_val):v(NULL),rows(r){
		v = allocator.allocate(r);
		std::uninitialized_fill_n(v, r, init_val);
	}
	TVector(const TVector &V):v(NULL),rows(V.size()){
		v = allocator.allocate(rows);
		std::uninitialized_copy(V.Raw(), V.Raw()+rows, v);
	}
	TVector(const TVectorView<value_type> &V):v(NULL),rows(V.size()){
		v = allocator.allocate(rows);
		std::uninitialized_copy(V.Raw(), V.Raw()+rows, v);
	}
	TVector& operator=(const TVector &V){
		Resize(V.size());
		for(size_t i = 0; i < rows; ++i){
			(*this)[i] = V[i];
		}
		return *this;
	}
	template <class VectorLikeType>
	TVector& operator=(const VectorLikeType &V){
		Resize(V.size());
		for(size_t i = 0; i < rows; ++i){
			(*this)[i] = V[i];
		}
		return *this;
	}
	virtual ~TVector(){
		if(NULL != v){ allocator.deallocate(v, rows); }
	}
	
	void Resize(size_t nRows){
		if(NULL != v){ allocator.deallocate(v, rows); }
		rows = nRows;
		v = allocator.allocate(rows);
	}
	
	size_t size() const{ return rows; }
	size_t Rows() const{ return size(); }
	NumericType operator[](size_t row) const{
#ifdef USE_MATRIX_OPS_ASSERTS
		assert(row < size());
#endif
		return v[row];
	}
	NumericType& operator[](size_t row){
#ifdef USE_MATRIX_OPS_ASSERTS
		assert(row < size());
#endif
		return v[row];
	}
	operator TVectorView<value_type>(){ return TVectorView<value_type>(v, rows, 1); }
	TVectorView<value_type> SubVector(size_t starting_row, size_t num_rows){
		return TVectorView<value_type>(&(v[starting_row]), num_rows, 1);
	}
	
	// Raw interface
	NumericType* Raw() const{ return v; }
	size_t Stride() const{ return 1; }
};

#endif // _TVECTOR_H_
