#ifndef _TVECTOR_H_
#define _TVECTOR_H_

#include "MatrixInterfaces.h"
#include <memory>
#include <cassert>
/*
// View for dense vectors
template <class T>
class TVectorView : public VectorViewBase<T>{
protected:
	T* V;
	size_t rows, stride;
public:
	typedef T value_type;
	typedef TVectorView<T> view_type;
	
	TVectorView(T* DataPtr, size_t nRows, size_t Stride):
		V(DataPtr),
		rows(nRows),
		stride(Stride){
	}
	virtual T& operator[](size_t row)      { return V[stride*row]; }
	virtual T  operator[](size_t row) const{ return V[stride*row]; }
	virtual size_t size() const{ return rows; }
};
*/

template <typename NumericType, class TAllocator = std::allocator<NumericType> >
class TVector : public WritableVector<NumericType>{
	NumericType *v;
	size_t rows;
	TAllocator allocator;
public:
	typedef NumericType value_type;
	typedef TVector<NumericType,TAllocator> self;
	
	typedef ReadableVector<value_type> non_view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVector<value_type> writable_vector;
	
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
		assert(row < size());
		return v[row];
	}
	NumericType& operator[](size_t row){
		assert(row < size());
		return v[row];
	}
	
	// Raw interface
	NumericType* Raw() const{ return v; }
	size_t Stride() const{ return 1; }
};


template <class T, class TAlloc>
class TrivialReadableVectorView<TVector<T,TAlloc> > : public ReadableVector<T>
{
	const TVector<T,TAlloc> &V;
public:
	typedef T value_type;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	
	TrivialReadableVectorView(const TVector<T,TAlloc> &vec):V(vec){}
	value_type operator[](size_t row) const{ return V[row]; }
	size_t size() const{ return V.size(); }
	
	value_type *Raw() const{ return V.Raw(); }
	size_t Stride() const{ return V.Stride(); }
};
template <class T, class TAlloc>
class TrivialWritableVectorView<TVector<T,TAlloc> > : public WritableVectorView<T>
{
	TVector<T,TAlloc> &V;
public:
	typedef T value_type;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVectorView<value_type> writable_vector;
	
	TrivialWritableVectorView(TVector<T,TAlloc> &vec):V(vec){}
	void Set(size_t row, const value_type &value) const{ V[row] = value; }
	const value_type& Get(size_t row) const{ return V[row]; }
	value_type& GetMutable(size_t row) const{ return V[row]; }
	size_t size() const{ return V.size(); }
	
	value_type *Raw() const{ return V.Raw(); }
	size_t Stride() const{ return V.Stride(); }
};






// Specialization of subvector view for dense vectors
template <typename T, class TAlloc>
class SubVectorView<TrivialWritableVectorView<TVector<T,TAlloc> > > : public WritableVectorView<T>{
	T* v;
	size_t rows, inc;
public:
	typedef T value_type;
	typedef TrivialWritableVectorView<TVector<T,TAlloc> > parent_view;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVectorView<value_type> writable_vector;
	
	SubVectorView<TrivialWritableVectorView<TVector<T,TAlloc> > >(TrivialWritableVectorView<TVector<T,TAlloc> > View, size_t RowStart, size_t nRows, size_t Inc = 1):v(&View.GetMutable(RowStart)),rows(nRows),inc(Inc*View.Stride()){}
	void Set(size_t row, const value_type &value) const{ v[inc*row] = value; }
	const value_type& Get(size_t row) const{ return v[inc*row]; }
	value_type& GetMutable(size_t row) const{ return v[inc*row]; }
	size_t size() const{ return rows; }
	
	size_t Stride() const{ return inc; }
	value_type *Raw() const{ return v; }
};
template <typename T, class TAlloc>
class SubVectorView<TrivialReadableVectorView<TVector<T,TAlloc> > > : public ReadableVector<T>{
	T* v;
	size_t rows, inc;
public:
	typedef T value_type;
	typedef TrivialReadableVectorView<TVector<T,TAlloc> > parent_view;
	typedef ReadableVector<value_type> view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVectorView<value_type> writable_vector;
	
	SubVectorView<TrivialReadableVectorView<TVector<T,TAlloc> > >(TrivialReadableVectorView<TVector<T,TAlloc> > View, size_t RowStart, size_t nRows, size_t Inc = 1):v(View.Raw()+RowStart),rows(nRows),inc(Inc*View.Stride()){}
	value_type operator[](size_t row) const{ return v[inc*row]; }
	size_t size() const{ return rows; }
	
	size_t Stride() const{ return inc; }
	value_type *Raw() const{ return v; }
};
template <class T, class TAlloc>
SubVectorView<TrivialWritableVectorView<TVector<T,TAlloc> > > SubVector(TVector<T> &V, size_t RowStart, size_t nRows, size_t inc = 1){
	return SubVectorView<TrivialWritableVectorView<TVector<T,TAlloc> > >(
		TrivialWritableVectorView<TVector<T,TAlloc> >(V),
		RowStart, nRows, inc);
}
template <class T, class TAlloc>
SubVectorView<TrivialReadableVectorView<TVector<T,TAlloc> > > SubVector(const TVector<T> &V, size_t RowStart, size_t nRows, size_t inc = 1){
	return SubVectorView<TrivialReadableVectorView<TVector<T,TAlloc> > >(
		TrivialReadableVectorView<TVector<T,TAlloc> >(V),
		RowStart, nRows, inc);
}

#endif // _TVECTOR_H_
