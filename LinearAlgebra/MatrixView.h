#ifndef _MATRIX_VIEW_H_
#define _MATRIX_VIEW_H_

template <class T>
class MatrixView{
protected:
	T* A;
	size_t rows, cols, col_stride;
public:
	MatrixView(T* DataPtr, size_t nRows, size_t nCols, size_t lda):
		A(DataPtr),
		rows(nRows),
		cols(nCols),
		col_stride(lda){
	}
	virtual T& operator()(size_t row, size_t col){
		return A[col*col_stride+row];
	}
	virtual T operator()(size_t row, size_t col) const{
		return A[col*col_stride+row];
	}
};

template <class T>
class MatrixViewView{
protected:
	MatrixView<T> view;
public:
	MatrixViewView(const MatrixView<T> &v):view(v){}
	virtual T& operator()(size_t row, size_t col){ return view(row,col); }
	virtual T operator()(size_t row, size_t col) const{ return view(row,col); }
};

template <class T>
class ConjugateTransposeView : public MatrixViewView<T>{
public:
	ConjugateTransposeView(const MatrixView<T> &view):MatrixViewView<T>(view){}
	T& operator()(size_t row, size_t col){ return std::conj(view(col, row)); }
	T operator()(size_t row, size_t col) const{ return std::conj(view(col, row)); }
};

template <class T>
class TransposeView : public MatrixViewView<T>{
public:
	TransposeView(const MatrixView<T> &view):MatrixViewView<T>(view){}
	T& operator()(size_t row, size_t col){ return view(col, row); }
	T operator()(size_t row, size_t col) const{ return view(col, row); }
};

#endif // _MATRIX_VIEW_H_
