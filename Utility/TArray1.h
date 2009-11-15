#ifndef _TARRAY1_H_
#define _TARRAY1_H_

#include <memory>
#include <cassert>

template <typename T, class TAllocator = std::allocator<T> >
class TArray1{
public:
	TArray1(size_t N, const T& init_val = T()):A(NULL),n(N){
		A = allocator.allocate(N);
		std::uninitialized_fill_n(A, N, init_val);
	}
	TArray1(const TArray1 &a):A(NULL),n(a.size()){
		A = allocator.allocate(n);
		std::uninitialized_copy(a.A, a.A+n, A);
	}
	TArray1& operator=(const TArray1 &a){
		if(this != &a){
			if(NULL != A){
				allocator.deallocate(A, n);
			}
			n = a.size();
			A = allocator.allocate(n);
			std::uninitialized_copy(a.A, a.A+n, A);
		}
		return *this;
	}
	~TArray1(){
		allocator.deallocate(A, n);
	}
	
	inline const T& operator[](size_t i) const{
		assert(i < n);
		return A[i];
	}
	inline T& operator[](size_t i){
		assert(i < n);
		return A[i];
	}
	inline size_t size() const{ return n; }
	inline T* Raw(){ return A; }
private:
	T* A;
	size_t n;
	TAllocator allocator;
};

#endif // _TARRAY1_H_
