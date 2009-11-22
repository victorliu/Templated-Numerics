#ifndef _MATRIX_INTERFACES_H_
#define _MATRIX_INTERFACES_H_

#include <cstddef> // for size_t

#include <cassert>

template <typename T>
class MatrixBase{
public:
	typedef T value_type;
};
template <typename T>
class VectorBase{
public:
	typedef T value_type;
};

// All actual matrices should implement either ReadableMatrix or
// WritableMatrix. All views should implement either ReadableMatrixView
// or WritableMatrixView. The same holds for vectors, correspondingly.

template <typename T>
class ReadableMatrix : public MatrixBase<T>{
public:
	typedef T value_type;
	virtual value_type operator()(size_t row, size_t col) const = 0;
	virtual value_type Get(size_t row, size_t col) const{ return (*this)(row,col); }
	virtual size_t Rows() const = 0;
	virtual size_t Cols() const = 0;
};

template <typename T>
class WritableMatrix : public ReadableMatrix<T>{
public:
	typedef T value_type;
	virtual value_type& operator()(size_t row, size_t col) = 0;
	virtual void Set(size_t row, size_t col, const value_type &value){ (*this)(row,col) = value; }
	value_type& GetMutable(size_t row, size_t col){ return (*this)(row,col); }
};

template <typename T>
class WritableMatrixView : public MatrixBase<T>{
public:
	typedef T value_type;
	value_type& operator()(size_t row, size_t col) const{ return GetMutable(row,col); }
	virtual void Set(size_t row, size_t col, const value_type &value) const = 0;
	virtual const value_type& Get(size_t row, size_t col) const = 0;
	virtual value_type& GetMutable(size_t row, size_t col) const = 0;
	virtual size_t Rows() const = 0;
	virtual size_t Cols() const = 0;
};








template <typename T>
class ReadableVector : public ReadableMatrix<T>, public VectorBase<T>{
public:
	typedef T value_type;
	virtual value_type operator[](size_t row) const = 0;
	virtual value_type Get(size_t row) const{ return (*this)[row]; }
	virtual size_t size() const = 0;
	
	size_t Rows() const{ return size(); }
	size_t Cols() const{ return 1; };
	
	value_type operator()(size_t row, size_t col) const{
		assert(0 == col);
		assert(row < Rows());
		return (*this)[row];
	}
};

template <typename T>
class WritableVector : public ReadableVector<T>{
public:
	typedef T value_type;
	virtual value_type& operator[](size_t row) = 0;
	virtual void Set(size_t row, const value_type &value){ (*this)[row] = value; }
	
	value_type& operator()(size_t row, size_t col){
		assert(0 == col);
		assert(row < ReadableVector<T>::Rows());
		return (*this)[row];
	}
};

template <typename T>
class WritableVectorView : public VectorBase<T>{
public:
	typedef T value_type;
	value_type& operator[](size_t row) const{ return GetMutable(row); }
	virtual void Set(size_t row, const value_type &value) const = 0;
	virtual const value_type& Get(size_t row) const = 0;
	virtual value_type& GetMutable(size_t row) const = 0;
	virtual size_t size() const = 0;
	
	size_t Rows() const{ return size(); }
	size_t Cols() const{ return 1; };
	
	void Set(size_t row, size_t col, const value_type &value) const{
		assert(0 == col);
		assert(row < ReadableVector<T>::Rows());
		return this->Set(row, value);
	}
	value_type Get(size_t row, size_t col) const{
		assert(0 == col);
		assert(row < ReadableVector<T>::Rows());
		return this->Get(row);
	}
};







// SFINAE selectors
template <class T, class R>
struct IsReadableMatrix{
};
template <class T, class R>
struct IsReadableMatrix<ReadableMatrix<T>, R>{
	typedef R type;
};

template <class T, class R>
struct IsWritableMatrix{
};
template <class T, class R>
struct IsWritableMatrix<WritableMatrix<T>, R>{
	typedef R type;
};

template <class T, class R>
struct IsNonWritableMatrix{
	typedef R type;
};
template <class T, class R>
struct IsNonWritableMatrix<WritableMatrix<T>, R>{
};
template <class T, class R>
struct IsNonWritableMatrix<WritableMatrixView<T>, R>{
};

template <class T, class R>
struct IsWritableMatrixView{
};
template <class T, class R>
struct IsWritableMatrixView<WritableMatrixView<T>, R>{
	typedef R type;
};



template <class T, class R>
struct IsReadableVector{
};
template <class T, class R>
struct IsReadableVector<ReadableVector<T>, R>{
	typedef R type;
};

template <class T, class R>
struct IsWritableVector{
};
template <class T, class R>
struct IsWritableVector<WritableVector<T>, R>{
	typedef R type;
};

template <class T, class R>
struct IsNonWritableVector{
	typedef R type;
};
template <class T, class R>
struct IsNonWritableVector<WritableVector<T>, R>{
};
template <class T, class R>
struct IsNonWritableVector<WritableVectorView<T>, R>{
};

template <class T, class R>
struct IsWritableVectorView{
};
template <class T, class R>
struct IsWritableVectorView<WritableVectorView<T>, R>{
	typedef R type;
};





template <class T, class R>
struct IsMatrixView{
};
template <class T, class R>
struct IsMatrixView<ReadableMatrix<T>, R>{
	typedef R type;
};
template <class T, class R>
struct IsMatrixNonView{
};
template <class T, class R>
struct IsMatrixNonView<ReadableMatrix<T>, R>{
	typedef R type;
};

template <class T, class R>
struct IsVectorView{
};
template <class T, class R>
struct IsVectorView<ReadableVector<T>, R>{
	typedef R type;
};
template <class T, class R>
struct IsVectorNonView{
};
template <class T, class R>
struct IsVectorNonView<ReadableVector<T>, R>{
	typedef R type;
};

#endif // _MATRIX_INTERFACES_H_
