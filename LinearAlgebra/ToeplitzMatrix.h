#ifndef _TOEPLITZ_MATRIX_H_
#define _TOEPLITZ_MATRIX_H_

#include <memory.h>
#include <vector>
/*
#include <iostream>
template <class value_type>
void pm(size_t N, const value_type *M, size_t ldm){
	std::cout << "{";
	for(int i = 0; i < N; ++i){
		std::cout << "{";
		std::cout << M[i+0*N];
		for(int j = 1; j < N; ++j){
			std::cout << ", " << M[i+j*N];
		}
		std::cout << "}," << std::endl;
	}
	std::cout << "}" << std::endl;
}
template <class value_type>
void pv(const std::vector<value_type> &v){
	std::cout << "{";
	std::cout << v[0];
	for(int i = 1; i < v.size(); ++i){
		std::cout << ", " << v[i];
	}
	std::cout << "}" << std::endl;
}
*/
template <class T>
class ToeplitzMatrix{
	typedef T value_type;
	// A toeplitz matrix (we store column major)
	// [ a0   a2n-2 a2n-3 ... an   ]
	// [ a1   a0    a2n-2     an+1 ]
	// [ a2   a1    a0        an+2 ]
	// [ ...  ...   ...       ...  ]
	// [ an-1 an-2  an-3  ... a0   ]
	size_t n;
	value_type *a; // length must be 2*n-1
public:
	ToeplitzMatrix(size_t rows):n(rows){
		a = new value_type[2*rows-1];
		memset(a, 0, (2*rows-1)*sizeof(value_type));
	}
	ToeplitzMatrix(size_t rows, value_type *fft):n(rows){
		a = new value_type[rows];
		memcpy(a, fft, (2*rows-1)*sizeof(value_type));
	}
	~ToeplitzMatrix(){
		delete [] a;
	}
	const value_type& operator[](size_t i) const{ return a[i]; }
	value_type& operator[](size_t i){ return a[i]; }
	const value_type& operator()(size_t row, size_t col) const{
		size_t N = 2*n-1;
		return a[(row-col+N)%N];
	}
	value_type& operator()(size_t row, size_t col){
		size_t N = 2*n-1;
		return a[(row-col+N)%N];
	}
	size_t size() const{ return n; }
	
	// Returns the inverse of this Toeplitz matrix in A, whose leading dimension is lda
	// Must have lda >= n
	// The inverse is stored in column-major format (Fortran style)
	int GetInverse(value_type *A, size_t lda) const{
		std::vector<value_type> rhop(n-1), rhon(n-1), e(n-1), g(n-1);
		value_type lambda;
		
		value_type ia0 = value_type(1)/a[0];
		for(size_t i = 1; i < n; ++i){
			rhop[i-1] = ia0*a[i];
			rhon[i-1] = ia0*a[2*n-1 - i]; // (2*n-1) - i
//std::cout << "rhop[" << i << "] = " << rhon[i] << std::endl;
		}
		
		// Initial values (iteration 1)
		lambda = 1 - rhon[0]*rhop[0];
		e[0] = -rhon[0];
		g[0] = -rhop[0];
//std::cout << "lambda = " << lambda << std::endl;
		for(size_t i = 1; i < n-1; ++i){
			// e and g are of length i at this point
			// Recursion
		
			// ahat is a reversed
			value_type e_dot_ahat(0);
			value_type r_dot_ghat(0);
			for(size_t j = 0; j < i; ++j){
//std::cout << "foo = " << e[j]<< std::endl;
//std::cout << "bar = " << rhon[n-2-j] << std::endl;
				e_dot_ahat += e[j]*rhon[i-1-j];
				r_dot_ghat += rhop[j]*g[i-1-j];
			}
//std::cout << "e.a^ = " << e_dot_ahat << std::endl;
//std::cout << "r.g^ = " << r_dot_ghat << std::endl;
			e[i] = -(rhon[i] + e_dot_ahat)/lambda;
			g[i] = -(rhop[i] + r_dot_ghat)/lambda;
//std::cout << "e[" << i << "] = " << e[i] << std::endl;
//std::cout << "g[" << i << "] = " << g[i] << std::endl;
			for(size_t j = 0; j < i; ++j){
				value_type ej = e[j];
				e[j] += e[i]*g[i-1-j];
				g[i-1-j] += g[i]*ej;
			}
			lambda -= e[i]*g[i]*lambda;
//std::cout << "e_" << i+1 << " = "; pv(e);
//std::cout << "g_" << i+1 << " = "; pv(g);
//std::cout << "lambda_" << i+1 << " = " << lambda << std::endl;
		}
		
		// Update the inverse
		// B is size (i+1) by (i+1) submatrix at lower right of A
		value_type ilambda(value_type(1)/lambda);
		value_type *B = A;//+(n-i-1)*lda+(n-i-1);
//std::cout << "inv lambda = " << ilambda << std::endl;
//std::cout << "e[0] = " << e[0] << std::endl;
		B[0+0*lda] = ilambda;
		for(size_t j = 1; j < n; ++j){
			B[0+j*lda] = e[j-1]*ilambda;
			B[j+0*lda] = g[j-1]*ilambda;
		}
//std::cout << "B = "; pm(i+1,B,lda);
		// Run over the triangle such that i+j < n, and i > 0 and j > 0
		for(size_t col = 1; col < n; ++col){
			for(size_t row = 1; row < n-col; ++row){
				B[row+col*lda] = B[row-1+(col-1)*lda] + ilambda*(g[row-1]*e[col-1]-e[(n-1)-row]*g[(n-1)-col]);
				// The persymmetric element to (row, col) is (n-col,n-row) for 0 based indexing
				if(row+col+1 != n){
					B[n-1-col+(n-1-row)*lda] = B[row+col*lda];
				}
			}
		}
		for(size_t col = 1; col < n; ++col){
			// copy the leading row and col
			B[n-1-col+(n-1-0)*lda] = B[0+col*lda];
			B[n-1-0+(n-1-col)*lda] = B[col+0*lda];
		}
		B[(n-1)+(n-1)*lda] = B[0];
//std::cout << "A = "; pm(n,A,lda);
	}
	void MultAdd(const value_type *x, value_type *y) const{ // y += this*x
	}
	// Returns the full Toeplitz matrix in A, whose leading dimension is lda
	// Must have lda >= n
	// The matrix is stored in column-major format (Fortran style)
	void Fill(value_type *A, size_t lda) const{ // A is filled in column-wise
		for(size_t i = 0; i < n; ++i){ // go through all columns of A
			size_t ni = n-i;
			memcpy(A+i*lda+i, a, ni*sizeof(value_type));
			memcpy(A+i*lda, a+n, i*sizeof(value_type));
		}
	}
};

#endif // _TOEPLITZ_MATRIX_H_
