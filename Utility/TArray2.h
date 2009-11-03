#ifndef _TARRAY2_H_
#define _TARRAY2_H_

#include <memory>

template <typename T, class TAllocator = std::allocator<T> >
class TArray2{
public:
	TArray2(size_t m, size_t n, const T& init_val = T()):A(NULL),n0(m),n1(n){
		const size_t mn = m*n;
		A = allocator.allocate(mn);
		std::uninitialized_fill_n(A, mn, init_val);
	}
	TArray2(const TArray2 &a):A(NULL),n0(a.Dim0()),n1(a.Dim1()){
		const size_t mn = n0*n1;
		A = allocator.allocate(mn);
		std::uninitialized_copy(a.A, a.A+mn, A);
	}
	TArray2& operator=(const TArray2 &a){
		if(this != &a){
			n0 = a.Dim0(); n1 = a.Dim1();
			const size_t mn = n0*n1;
			A = allocator.allocate(mn);
			std::uninitialized_copy(a.A, a.A+mn, A);
		}
		return *this;
	}
	~TArray2(){
		allocator.deallocate(A, n0*n1);
	}
	
	inline const T& operator()(size_t i, size_t j) const{ return A[i*n1+j]; }
	inline       T& operator()(size_t i, size_t j)      { return A[i*n1+j]; }
	inline size_t Dim0() const{ return n0; }
	inline size_t Dim1() const{ return n1; }
	inline T* Raw(){ return A; }
private:
	T* A;
	size_t n0, n1;
	TAllocator allocator;
};

#endif // _TARRAY2_H_
