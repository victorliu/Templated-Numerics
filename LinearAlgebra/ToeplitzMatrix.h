#ifndef _TOEPLITZ_MATRIX_H_
#define _TOEPLITZ_MATRIX_H_

#include <memory.h>
#include <vector>

#ifdef USE_MATRIX_VIEW
# include "TMatrix.h"
#endif // USE_MATRIX_VIEW


// Optional preprocessor switches:
//   USE_TMATRIX_BASE_INTERFACE
//   USE_MATRIX_ASSERTS
//   USE_MATRIX_VIEW

#ifdef USE_MATRIX_ASSERTS
# include <cassert>
#endif
#include <iostream>
/*template <class value_type>
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
class ToeplitzMatrix
#ifdef USE_TMATRIX_BASE_INTERFACE
	: public TMatrixBase<T>
#endif
{
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
	// For an FFT data set of odd length, we don'tneed to pad
	// For an even length FFT, we need to pad an extra element
	ToeplitzMatrix(size_t fft_len, value_type *fft, size_t fft_stride = 1):n(fft_len/2+1){
		a = new value_type[2*n-1];
		if(1 == fft_stride){
			memcpy(a, fft, fft_len*sizeof(value_type));
		}else{
			for(size_t i = 0; i < fft_len; ++i){
				a[i] = fft[fft_stride*i];
			}
		}
		if(0 == (fft_len & 1)){
			a[2*n-2] = value_type(0);
		}
	}
	~ToeplitzMatrix(){
		delete [] a;
	}
	const value_type& operator[](size_t i) const{ return a[i]; }
	value_type& operator[](size_t i){ return a[i]; }
	value_type operator()(size_t row, size_t col) const{
		size_t N = 2*n-1;
		return a[(row-col+N)%N];
	}
	value_type& operator()(size_t row, size_t col){
		size_t N = 2*n-1;
		return a[(row-col+N)%N];
	}
	size_t size() const{ return n; }
	size_t Rows() const{ return n; }
	size_t Cols() const{ return n; }
	
	class InverseMatrix
#ifdef USE_TMATRIX_BASE_INTERFACE
		: public TMatrixBase<T>
#endif
	{
		size_t n;
		value_type *lambda;
		value_type *g, *e;
	public:
		InverseMatrix(const ToeplitzMatrix &A):n(A.size()){
			lambda = new value_type[2*(n-1)+1];
			g = lambda + 1;
			e = g+(n-1);
			
			std::vector<value_type> rhop(n-1), rhon(n-1);
			
			value_type ia0 = value_type(1)/A(0,0);
			for(size_t i = 1; i < n; ++i){
				rhop[i-1] = ia0*A(i,0);
				rhon[i-1] = ia0*A(0,i); // (2*n-1) - i
	//std::cout << "rhop[" << i << "] = " << rhon[i] << std::endl;
			}
			
			// Initial values (iteration 1)
			*lambda = value_type(1) - rhon[0]*rhop[0];
			e[0] = -rhon[0];
			g[0] = -rhop[0];
	//std::cout << "lambda = " << *lambda << std::endl;
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
				e[i] = -(rhon[i] + e_dot_ahat)/(*lambda);
				g[i] = -(rhop[i] + r_dot_ghat)/(*lambda);
	//std::cout << "e[" << i << "] = " << e[i] << std::endl;
	//std::cout << "g[" << i << "] = " << g[i] << std::endl;
				for(size_t j = 0; j < i; ++j){
					value_type ej = e[j];
					e[j] += e[i]*g[i-1-j];
					g[i-1-j] += g[i]*ej;
				}
				(*lambda) -= e[i]*g[i]*(*lambda);
	//std::cout << "e_" << i+1 << " = "; pv(e);
	//std::cout << "g_" << i+1 << " = "; pv(g);
	//std::cout << "lambda_" << i+1 << " = " << lambda << std::endl;
			}
			
			// Update the inverse
			// B is size (i+1) by (i+1) submatrix at lower right of A
			*lambda /= ia0; // put the factor of a0 back in
			value_type ilambda = value_type(1)/(*lambda);
			for(size_t i = 0; i < n-1; ++i){
				g[i] *= ilambda;
				e[i] *= ilambda;
			}
		}
		~InverseMatrix(){
			delete [] lambda;
		}
		size_t Rows() const{ return n; }
		size_t Cols() const{ return n; }
		void Fill(value_type *A, size_t lda){
	//std::cout << "inv lambda = " << ilambda << std::endl;
	//std::cout << "e[0] = " << e[0] << std::endl;
			A[0+0*lda] = value_type(1)/(*lambda);
			for(size_t j = 1; j < n; ++j){
				A[0+j*lda] = e[j-1];
				A[j+0*lda] = g[j-1];
			}
	//std::cout << "A = "; pm(i+1,A,lda);
			// Run over the triangle such that i+j < n, and i > 0 and j > 0
			for(size_t col = 1; col < n; ++col){
				for(size_t row = 1; row < n-col; ++row){
					A[row+col*lda] = A[row-1+(col-1)*lda] + (*lambda)*(g[row-1]*e[col-1]-e[(n-1)-row]*g[(n-1)-col]);
					// The persymmetric element to (row, col) is (n-col,n-row) for 0 based indexing
					if(row+col+1 != n){
						A[n-1-col+(n-1-row)*lda] = A[row+col*lda];
					}
				}
			}
			for(size_t col = 1; col < n; ++col){
				// copy the leading row and col
				A[n-1-col+(n-1-0)*lda] = A[0+col*lda];
				A[n-1-0+(n-1-col)*lda] = A[col+0*lda];
			}
			A[(n-1)+(n-1)*lda] = A[0];
	//std::cout << "A = "; pm(n,A,lda);
		}
		size_t size() const{ return n; }
		value_type& operator()(size_t Row, size_t Col){
			static value_type black_hole;
#ifdef USE_MATRIX_ASSERTS
			assert(0);
#endif
			return black_hole;
		}
		value_type operator()(size_t Row, size_t Col) const{ // this is broken
			if(Col > n-1-Row){
				size_t t = n-1-Col;
				Col = n-1-Row;
				Row = t;
			}
			
			std::vector<value_type> Brow(n);
			// Fill in the first row
			Brow[0] = value_type(1)/(*lambda);
			for(size_t col = 1; col <= Col; ++col){
				Brow[col] =  e[col-1];
			}
			for(size_t row = 1; row <= Row; ++row){
				// Compute the new row
				value_type prev = Brow[0];
				for(size_t col = 1; col <= Col; ++col){
					value_type cur = prev + (*lambda)*(g[row-1]*e[col-1]-e[(n-1)-row]*g[(n-1)-col]);
					prev = Brow[col];
					Brow[col] = cur;
				}
				Brow[0] = g[row-1];
			}
			return Brow[Col];
		}
		// C = alpha*this*B, B is size n by M
		void MultAdd(size_t M, const value_type *B, size_t ldb,
			value_type *C, size_t ldc,
			const value_type &alpha = value_type(1)) const{
			
			for(size_t row = 0; row < n; ++row){
				for(size_t col = 0; col < M; ++col){
					value_type sum(0);
					for(size_t mid = 0; mid < n; ++mid){
						sum += ((*this)(row,mid)) * B[mid+ldb*col];
					}
					C[row+ldc*col] = sum;
				}
			}
			return;
			
			// First row
			std::vector<value_type> Brow(n);
			// Fill in the first row
			Brow[0] = alpha/(*lambda);
			for(size_t i = 0; i < n-1; ++i){
				Brow[i+1] =  alpha * e[i];
			}
			// Update using the first/last element
			{size_t row = 0;
//for(size_t i = 0; i < n-row; ++i){ std::cout << Brow[i] << " "; }std::cout << std::endl;
				{size_t mid = 0;
					for(size_t col = 0; col < M; ++col){
						C[row+ldc*col] += Brow[mid] * B[mid+ldb*col];
						C[(n-1-mid)+ldc*col] += Brow[mid] * B[(n-1-row)+ldb*col];
					}
				}
			}
			// Update using the first row/last column
			{size_t row = 0;
				for(size_t mid = 1; mid < n-1; ++mid){
					for(size_t col = 0; col < M; ++col){
						C[row+ldc*col] += Brow[mid] * B[mid+ldb*col];
						C[(n-1-mid)+ldc*col] += Brow[mid] * B[(n-1-row)+ldb*col];
					}
				}
				// Update the anti-diagonal
				{size_t mid = n-1-row;
					for(size_t col = 0; col < M; ++col){
						C[row+ldc*col] += Brow[mid] * B[mid+ldb*col];
					}
				}
			}
			for(size_t row = 1; row < n; ++row){
				// Compute the new row
				value_type prev = Brow[0];
				for(size_t i = 1; i <= n-row; ++i){
					value_type cur = prev + alpha * (*lambda)*(g[row-1]*e[i-1]-e[(n-1)-row]*g[(n-1)-i]);
					prev = Brow[i];
					Brow[i] = cur;
				}
				Brow[0] = alpha * g[row-1];
				
//for(size_t i = 0; i < n-row; ++i){ std::cout << Brow[i] << " "; }std::cout << std::endl;
				// Update using the first column/last row
				if(row < n-1){ size_t mid = 0;
					for(size_t col = 0; col < M; ++col){
						C[row+ldc*col] += Brow[mid] * B[mid+ldb*col];
						C[(n-1-mid)+ldc*col] += Brow[mid] * B[(n-1-row)+ldb*col];
					}
				}
				// Update the interior
				for(size_t mid = 1; mid < n-1-row; ++mid){ // go up to just before the anti-diagonal
					for(size_t col = 0; col < M; ++col){
						C[row+ldc*col] += Brow[mid] * B[mid+ldb*col];
						C[(n-1-mid)+ldc*col] += Brow[mid] * B[(n-1-row)+ldb*col];
					}
				}
				// Update the anti-diagonal
				{size_t mid = n-1-row;
					for(size_t col = 0; col < M; ++col){
						C[row+ldc*col] += Brow[mid] * B[mid+ldb*col];
					}
				}
			}
		}
	};
	// Returns the inverse of this Toeplitz matrix in A, whose leading dimension is lda
	// Must have lda >= n
	// The inverse is stored in column-major format (Fortran style)
	int GetInverse(value_type *A, size_t lda) const{
		InverseMatrix i(*this);
		i.Fill(A, lda);
		return 0;
	}
	// C = alpha*this*B + beta*C, B is size n by M
	void MultAdd(size_t M, const value_type *B, size_t ldb,
		value_type *C, size_t ldc,
		const value_type &alpha = value_type(1),
		const value_type &beta = value_type(0)) const{
		for(size_t row = 0; row < n; ++row){
			for(size_t col = 0; col < M; ++col){
				value_type sum(0);
				for(size_t mid = 0; mid < n; ++mid){
					sum += ((*this)(row,mid)) * B[mid+ldb*col];
				}
				C[row+ldc*col] = alpha * sum + beta * C[row+ldc*col];
			}
		}
	}
	
	// Returns the full Toeplitz matrix in A, whose leading dimension is lda
	// Must have lda >= n
	// The matrix is stored in column-major format (Fortran style)
	void Fill(value_type *A, size_t lda) const{ // A is filled in column-wise
		for(size_t i = 0; i < n; ++i){ // go through all columns of A
			size_t ni = n-i;
			memcpy(A+i*lda+i, a, ni*sizeof(value_type)); // copy lower triangular part
			for(size_t j = i+1; j < n; ++j){
				A[i+lda*j] = a[2*n-1-(j-i)];
			}
		}
	}
	
#ifdef USE_MATRIX_VIEW
public:
	class ToeplitzMatrixView : public MatrixViewBase<value_type>{
		value_type *a;
		size_t n;
	public:
		typedef MatrixViewBase<value_type> value_type;
		typedef ToeplitzMatrixView<value_type> matrix_type
		
		ToeplitzMatrixView(value_type *DataPtr, size_t size):a(DataPtr),n(size){}
		value_type& operator()(size_t row, size_t col){
			size_t N = 2*n-1;
			return a[(row-col+N)%N];
		}
		value_type operator()(size_t row, size_t col) const{
			size_t N = 2*n-1;
			return a[(row-col+N)%N];
		}
		size_t Rows() const{ return n; }
		size_t Cols() const{ return n; }
	};
#endif // USE_MATRIX_VIEW
};

#endif // _TOEPLITZ_MATRIX_H_
