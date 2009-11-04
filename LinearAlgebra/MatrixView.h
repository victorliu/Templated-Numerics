#ifndef _MATRIX_VIEW_H_
#define _MATRIX_VIEW_H_

template <class T>
class VectorViewBase{
public:
	typedef T value_type;
	//typedef [self]<T> vector_type
	//typedef ... parent_view;
	virtual value_type  operator[](size_t row) const = 0;
	virtual value_type& operator[](size_t row)       = 0;
	virtual size_t size() const = 0;
	size_t Rows() const{ return size(); }
	size_t Cols() const{ return 1; }
};

template <class ViewClass>
class ConjugateView : public VectorViewBase<typename ViewClass::value_type>{
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	ConjugateView(const ViewClass &view):ViewClass(view){}
	value_type  operator[](size_t row) const{ return std::conj(view[row]); }
	value_type& operator[](size_t row)      { return std::conj(view[row]); } // illegal to use
	size_t size() const{ return view.size(); }
};

template <class ViewClass>
class SubVectorView : public VectorViewBase<typename ViewClass::value_type>{
	ViewClass view;
	size_t row_start, rows;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	SubVectorView(const ViewClass &view, size_t RowStart, size_t nRows):ViewClass(view),row_start(RowStart),rows(nRows){}
	value_type  operator[](size_t row) const{ return view[row-row_start]; }
	value_type& operator[](size_t row)      { return view[row-row_start]; }
	size_t size() const{ return rows; }
};


template <class T>
class MatrixViewBase{
public:
	typedef T value_type;
	//typedef ... parent_view;
	//typedef [self]<T> matrix_type
	virtual value_type  operator()(size_t row, size_t col) const = 0;
	virtual value_type& operator()(size_t row, size_t col)       = 0;
	virtual size_t Rows() const = 0;
	virtual size_t Cols() const = 0;
};

template <class ViewClass>
class TransposeView : public MatrixViewBase<typename ViewClass::value_type>{
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	TransposeView(const ViewClass &view):ViewClass(view){}
	value_type  operator()(size_t row, size_t col) const{ return view(col, row); }
	value_type& operator()(size_t row, size_t col)      { return view(col, row); }
	size_t Rows() const{ return view.Cols(); }
	size_t Cols() const{ return view.Rows(); }
};

template <class ViewClass>
class ConjugateTransposeView : public MatrixViewBase<typename ViewClass::value_type>{
	ViewClass view;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	ConjugateTransposeView(const ViewClass &view):ViewClass(view){}
	value_type  operator()(size_t row, size_t col) const{ return std::conj(view(col, row)); }
	value_type& operator()(size_t row, size_t col)      { return std::conj(view(col, row)); } // illegal to use
	size_t Rows() const{ return view.Cols(); }
	size_t Cols() const{ return view.Rows(); }
};

template <class ViewClass>
class SubMatrixView : public MatrixViewBase<typename ViewClass::value_type>{
	ViewClass view;
	size_t row_start, col_start, rows, cols;
public:
	typedef typename ViewClass::value_type value_type;
	typedef ViewClass parent_view;
	SubMatrixView(const ViewClass &view, size_t RowStart, size_t ColStart, size_t nRows, size_t nCols):ViewClass(view),row_start(RowStart),col_start(ColStart),rows(nRows),cols(nCols){}
	value_type  operator()(size_t row, size_t col) const{ return view(row-row_start, col-col_start); }
	value_type& operator()(size_t row, size_t col)      { return view(row-row_start, col-col_start); }
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
};

#endif // _MATRIX_VIEW_H_
