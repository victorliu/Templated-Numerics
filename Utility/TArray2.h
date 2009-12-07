#ifndef _TARRAY2_H_
#define _TARRAY2_H_

#include <memory>

template <typename T, class TAllocator = std::allocator<T> >
class TArray2{
public:
	TArray2():A(NULL),n0(0),n1(0){}
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
			if(NULL != A){
				allocator.deallocate(A, n0*n1);
			}
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
	
	void Resize(size_t N0, size_t N1){
		if(N0 == n0 && N1 == n1){ return; }
		allocator.deallocate(A, n0*n1);
		n0 = N0; n1 = N1;
		A = allocator.allocate(n0*n1);
	}
	void Resize(size_t N0, size_t N1, const T &fill_value){
		if(N0 == n0 && N1 == n1){ return; }
		allocator.deallocate(A, n0*n1);
		n0 = N0; n1 = N1;
		A = allocator.allocate(n0*n1);
		std::uninitialized_fill_n(A, n0*n1, fill_value);
	}

	// Internally we use the C-style multidimensional array indexing
	// where the last index changes fastest.
	inline size_t Idx(size_t i, size_t j) const{ return i*n1+j; }
	inline void UnIdx(size_t idx, size_t &i, size_t &j){ i = idx/n1; j = idx%n1; }
	inline const T& operator()(size_t i, size_t j) const{ return A[Idx(i,j)]; }
	inline       T& operator()(size_t i, size_t j)      { return A[Idx(i,j)]; }
	inline size_t Dim0() const{ return n0; }
	inline size_t Dim1() const{ return n1; }
	inline T* Raw(){ return A; }
private:
	T* A;
	size_t n0, n1;
	TAllocator allocator;
};

#endif // _TARRAY2_H_
