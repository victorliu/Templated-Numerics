#include <iostream>
#define USE_MATRIX_OPS_ASSERTS
#define USE_COMPLEX_MATRICES
#define USE_ADVANCED_MATRIX_OPS
#include "TMatrix.h"
#include "TVector.h"
#include "TPermutationMatrix.h"
#include "MatrixViews.h"
#include "MatrixIO.h"
#include "MatrixOps.h"
#include "../NumericTypes/complex_io.hpp"

using namespace MatrixOps;

typedef double real_t;
typedef std::complex<real_t> complex_t;

typedef TMatrix<complex_t> matrix_t;
typedef TVector<complex_t> vector_t;
typedef TPermutationMatrix<size_t> perm_t;

int main(){
	matrix_t M(3,3), A(3,3), C(3,3);
	vector_t x(3);
	
	
	Fill(x, 0);
	Fill(M, 0);
	
	int count = 0;
	for(size_t j = 0; j < 3; ++j){
		for(size_t i = 0; i < 3; ++i){
			M(i,j) = matrix_t::value_type(count++,count++);
		}
	}
	for(size_t i = 0; i < 3; ++i){
		x[i] = vector_t::value_type(count++,count++);
	}
	std::cout << M << std::endl;
	std::cout << x << std::endl;
	
	
	vector_t y(3);
	Mult(ConjugateTranspose(Transpose(M)), x, y);
	std::cout << y << std::endl;
	
	std::cout << Diagonal(M) << std::endl;
	
	Copy(x, y);
	Copy(M, A);
	
	std::cout << Dot(x,y) << std::endl;
	std::cout << ConjugateDot(x,y) << std::endl;
	
	Scale(x, matrix_t::value_type(2));
	Scale(A, matrix_t::value_type(3));
	
	Add(x, y, matrix_t::value_type(2));
	Add(A, M, matrix_t::value_type(2));
	
	Rank1Update(A, x, y);
	Rank1Update(A, x);
	Rank2Update(A, x, y);
	Mult(A, M, C);
	
	std::cout << LargestElementIndex(x) << std::endl;

	A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
	A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;
	A(2,0) = 3; A(2,1) = 8; A(2,2) = 9;

	y[0] = 1; y[1] = 1; y[2] = 1;
	Solve(A, y, x);
	std::cout << A << std::endl;
	std::cout << x << std::endl;
	std::cout << y << std::endl;
	
	return 0;
}
