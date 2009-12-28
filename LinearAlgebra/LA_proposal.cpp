#include <iostream>

template <typename T>
class MatrixBase{
public:
	typedef T value_type;
};

template <typename T>
class WritableMatrix : public MatrixBase<T>{
public:
	typedef T value_type;
	virtual const value_type& operator()() const = 0;
	virtual value_type& operator()() = 0;
};

template <typename T>
class NonwritableMatrix : public MatrixBase<T>{
public:
	typedef T value_type;
	virtual value_type operator()() const = 0;
};



//struct IsRectangularMatrix
//struct IsVector
//struct IsTriangularMatrix

template <class T, class R>
struct IsWritableMatrix{
};
template <class T, class R>
struct IsWritableMatrix<WritableMatrix<T>, R>{
	typedef R type;
};
template <class T, class R>
struct IsNonwritableMatrix{
};
template <class T, class R>
struct IsNonwritableMatrix<NonwritableMatrix<T>, R>{
	typedef R type;
};

template <class Parent>
class TrivialView : public WritableMatrix<typename Parent::value_type>{
	Parent& parent;
public:
	typedef typename Parent::value_type value_type;
	typedef WritableMatrix<value_type> writability;
	
	TrivialView(const Parent &parent):parent(const_cast<Parent&>(parent)){}
	const value_type& operator()() const{ return parent(); }
	value_type& operator()(){
		return parent();
	}
};



template <class T>
class TMatrix : public WritableMatrix<T>{
	T m;
public:
	typedef T value_type;
	typedef WritableMatrix<T> writability;
	
	TMatrix():m(0){}
	const value_type& operator()() const{ return m; }
	value_type& operator()(){ return m; }
	
	operator TrivialView<TMatrix<T> >(){
		return TrivialView<TMatrix<T> >(*this);
	}
};





template <bool B, class TrueType, class FalseType>
struct IfThenElse{
};
template <class TrueType, class FalseType>
struct IfThenElse<true, TrueType, FalseType>{
	typedef TrueType type;
};
template <class TrueType, class FalseType>
struct IfThenElse<false, TrueType, FalseType>{
	typedef FalseType type;
};



template <class Parent, bool RO>
class TransposeView : public IfThenElse<RO, NonwritableMatrix<typename Parent::value_type>, WritableMatrix<typename Parent::value_type> >::type{
	typename IfThenElse<RO, const Parent, Parent>::type *parent;
public:
	typedef typename Parent::value_type value_type;
	typedef typename IfThenElse< RO, NonwritableMatrix<typename Parent::value_type>, WritableMatrix<typename Parent::value_type> >::type writability;
	
	TransposeView(const Parent &parent):parent(const_cast<typename IfThenElse<RO, const Parent, Parent>::type *>(&parent)){}
	
	value_type& operator()(){ return (*parent)(); }
	
	typename IfThenElse<RO,value_type,const value_type&>::type operator()() const{ return (*parent)(); }
};

template <class Parent>
class ScaledView : public NonwritableMatrix<typename Parent::value_type>{
	const Parent &parent;
	const typename Parent::value_type scale;
public:
	typedef typename Parent::value_type value_type;
	typedef NonwritableMatrix<typename Parent::value_type> writability;
	
	ScaledView(const Parent &parent, const value_type &scale):parent(parent),scale(scale){}
	
	value_type operator()() const{ return scale*parent(); }
};




template <class T, class R>
struct IsConst{
};
template <class T, class R>
struct IsConst<const T, R>{
	typedef R type;
};
template <class T, class R>
struct IsNonConst{
	typedef R type;
};
template <class T, class R>
struct IsNonConst<const T, R>{
};




template <class M>
typename IsConst<M,
TransposeView<M,true> >::type Transpose(const M &parent){
	return TransposeView<M,true>(parent);
}
template <class M>
typename IsNonwritableMatrix<typename M::writability,
typename IsNonConst<M,
TransposeView<M,true>
	>::type>::type
Transpose(const M &parent){
	return TransposeView<M,true>(parent);
}
template <class M>
typename IsWritableMatrix<typename M::writability,
typename IsNonConst<M,
TransposeView<M,false>
	>::type>::type
Transpose(const M &parent){
	return TransposeView<M,false>(parent);
}



template <class M>
ScaledView<M> Scale(const M &parent, const typename M::value_type &scale){
	return ScaledView<M>(parent, scale);
}




template <class SRC, class DST>
typename IsWritableMatrix<typename DST::writability,
int
	>::type
Copy(const SRC &src, const DST &dst_){
	DST &dst = const_cast<DST&>(dst_);
	dst() = src();
	return 0;
}

template <bool SrcRO>
void Copy(const TransposeView<typename IfThenElse<SrcRO,const TMatrix<int>, TMatrix<int> >::type,SrcRO> &src, TMatrix<int> &dst){
	std::cout << "In specialization const" << std::endl;
	dst() = src();
}


int main(){
	const TMatrix<int> A;
	TMatrix<int> B;
	
	B() = 1;
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	Copy(A, Transpose(B));
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	
	
	B() = 1;
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	Copy(Transpose(A), B);
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	
	
	B() = 1;
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	Copy(A, B);
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	
	
	B() = 1;
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	Copy(Transpose(A), Transpose(B));
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	
	
	B() = 1;
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	Copy(Transpose(A), Transpose(Transpose(B)));
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	
	B() = 1;
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	Copy(Transpose(Transpose(A)), Transpose(B));
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	
	
	B() = 1;
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	Copy(Transpose(A), Transpose(B));
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	
	
	B() = 1;
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	Copy(Transpose(Transpose(A)), Transpose(Transpose(B)));
	std::cout << "A: " << A() << ", B: " << B() << std::endl;
	
	
	TMatrix<int> C;
	TMatrix<int> D;
	
	D() = 1;
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	Copy(C, Transpose(D));
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	
	
	D() = 1;
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	Copy(Transpose(C), D);
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	
	
	D() = 1;
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	Copy(C, D);
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	
	
	D() = 1;
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	Copy(Transpose(C), Transpose(D));
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	
	
	
	
	
	
	
	D() = 1;
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	Copy(Transpose(C), Transpose(Transpose(D)));
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	
	
	D() = 1;
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	Copy(Transpose(Transpose(C)), Transpose(D));
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	
	
	D() = 1;
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	Copy(Transpose(C), Transpose(D));
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	
	
	D() = 1;
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	Copy(Transpose(Transpose(C)), Transpose(Transpose(D)));
	std::cout << "C: " << C() << ", D: " << D() << std::endl;
	
	return 0;
}



