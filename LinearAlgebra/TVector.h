#ifndef _TVECTOR_H_
#define _TVECTOR_H_

#include "MatrixView.h"

template <typename NumericType>
class TVectorBase{
public:
	typedef NumericType value_type;
	// Implementers must define this:
	//typedef VectorViewBase<T> View;
	
	virtual size_t size() const = 0;
	size_t Rows() const{ return size(); }
	size_t Cols() const{ return 1; }
	virtual const NumericType& operator[](size_t row) const = 0;
	virtual       NumericType& operator[](size_t row)       = 0;
	
	//virtual operator VectorViewBase&() = 0;
};



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
	
	friend class SubVectorView<TVectorView>;
};

// Specialization of submatrix view for dense matrices
template <typename T>
class SubVectorView<TVectorView<T> > : public VectorViewBase<typename TVectorView<T>::value_type>{
	TVectorView view;
public:
	typedef typename TVectorView<T>::value_Type value_type;
	SubVectorView(const TVectorView<T> &view, size_t RowStart, size_t nRows):TVectorView<T>(&view[RowStart], nRows, view.stride){}
	value_type  operator[](size_t row) const{ return view[row]; }
	value_type& operator[](size_t row)      { return view[row]; }
	size_t size() const{ return view.size(); }
};


template <typename NumericType, class TAllocator = std::allocator<NumericType> >
class TVector : public TVectorBase<NumericType>{
	NumericType *v;
	size_t rows;
	TAllocator allocator;
public:
	typedef NumericType value_type;
	typedef TVectorView<NumericType> View;
	
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
	TVector& operator=(const TVector &V){
		if(this != &V){
			Resize(V.size()); // cannot update rows and cols yet
			std::uninitialized_copy(V.Raw(), V.Raw()+rows, v);
		}
		return *this;
	}
	TVector(const TVectorView<value_type> &V):v(NULL),rows(V.size()){
		v = allocator.allocate(rows);
		std::uninitialized_copy(V.Raw(), V.Raw()+rows, v);
	}
	template <class VectorLikeType>
	TVector(const VectorLikeType &V):v(NULL),rows(V.size()){
		v = allocator.allocate(rows);
		for(size_t i = 0; i < rows; ++i){
			(*this)[i] = V[i];
		}
	}
	~TVector(){
		allocator.deallocate(v, rows);
	}
	
	void Resize(size_t nRows){
		allocator.deallocate(v, rows);
		rows = nRows;
		v = allocator.allocate(rows);
	}
	
	size_t size() const{ return rows; }
	size_t Rows() const{ return size(); }
	const NumericType& operator[](size_t row) const{ return v[row]; }
	NumericType& operator[](size_t row){ return v[row]; }
	operator TVectorView<value_type>(){ return TVectorView<value_type>(v, rows, 1); }
	TVectorView<value_type> SubVector(size_t starting_row, size_t num_rows){
		return TVectorView<value_type>(&(v[starting_row]), num_rows, 1);
	}
	
	// Raw interface
	NumericType* Raw() const{ return v; }
	size_t Stride() const{ return 1; }
};

#endif // _TVECTOR_H_
