#ifndef _MATRIX_VIEW_H_
#define _MATRIX_VIEW_H_

template <class T>
class MatrixViewBase{
protected:
	T* A;
	size_t rows, cols, col_stride;
public:
	MatrixViewBase(T* DataPtr, size_t nRows, size_t nCols, size_t lda):
		A(DataPtr),
		rows(nRows),
		cols(nCols),
		col_stride(lda){
	}
	T& operator()(size_t row, size_t col){
		return A[col*col_stride+row];
	}
	const T& operator()(size_t row, size_t col) const{
		return A[col*col_stride+row];
	}
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
};

template <class T>
class MatrixView{
protected:
	MatrixViewBase<T> view;
public:
	MatrixView(const MatrixViewBase<T> &v):view(v){}
	virtual T& operator()(size_t row, size_t col){ return view(row,col); }
	virtual T operator()(size_t row, size_t col) const{ return view(row,col); }
	virtual size_t Rows() const{ return view.Rows(); }
	virtual size_t Cols() const{ return view.Cols(); }
};

template <class T>
class ConjugateTransposeView : public MatrixView<T>{
	MatrixViewBase<T> view;
public:
	ConjugateTransposeView(const MatrixViewBase<T> &view):MatrixView<T>(view){}
	T& operator()(size_t row, size_t col){ return conj(view(col, row)); }
	T operator()(size_t row, size_t col) const{ return conj(view(col, row)); }
	size_t Rows() const{ return view.Cols(); }
	size_t Cols() const{ return view.Rows(); }
};

template <class T>
class TransposeView : public MatrixView<T>{
	MatrixViewBase<T> view;
public:
	TransposeView(const MatrixViewBase<T> &view):MatrixView<T>(view){}
	T& operator()(size_t row, size_t col){ return view(col, row); }
	T operator()(size_t row, size_t col) const{ return view(col, row); }
	size_t Rows() const{ return view.Cols(); }
	size_t Cols() const{ return view.Rows(); }
};

#endif // _MATRIX_VIEW_H_
