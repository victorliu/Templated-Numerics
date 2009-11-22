#include <iostream>
#define USE_COMPLEX_MATRICES
#define USE_ADVANCED_MATRIX_OPS
#include "TMatrix.h"
#include "TIdentityMatrix.h"
#include "TVector.h"
#include "MatrixIO.h"
#include "MatrixOps.h"
#include "../NumericTypes/complex_io.hpp"

using namespace MatrixOps;

typedef double real_t;
typedef std::complex<real_t> complex_t;

typedef TMatrix<complex_t> matrix_t;
typedef TIdentityMatrix<complex_t> identity_matrix_t;
typedef TVector<complex_t> vector_t;

int main(){
	matrix_t M(3,3), A(3,3), C(3,3);
	vector_t x(3), y(3), z(3);
	
	Copy(identity_matrix_t(3), A);
	std::cout << A << std::endl;
	
	Copy(Scaled(identity_matrix_t(3), complex_t(2)), A);
	std::cout << A << std::endl;
	Copy(Scaled(identity_matrix_t(2), complex_t(4)), SubMatrix(A,1,1,2,2));
	std::cout << A << std::endl;
	Copy(Scaled(Scaled(identity_matrix_t(2), complex_t(2)), complex_t(16)), SubMatrix(A,1,1,2,2));
	std::cout << A << std::endl;
	Copy(Scaled(Scaled(identity_matrix_t(1), complex_t(4)), complex_t(2)), SubMatrix(SubMatrix(A,1,1,2,2), 0,0,1,1));
	std::cout << A << std::endl;
	
	Copy(Diagonal(A), x);
	std::cout << x << std::endl;
	Copy(x, y);
	std::cout << y << std::endl;
	Copy(SubMatrix(A,1,1,2,2), SubMatrix(C,0,0,2,2));
	std::cout << A << std::endl;
	Copy(y, Diagonal(C));
	std::cout << C << std::endl;
	
	Fill(Diagonal(SubMatrix(A,0,0,3,3)), complex_t(1));
	std::cout << A << std::endl;
	
	Copy(x, y);
	Fill(x, complex_t(1));
	
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
	
	
	Mult(ConjugateTranspose(Transpose(M)), x, y);
	std::cout << y << std::endl;
	
	std::cout << Diagonal(M) << std::endl;
	
	Copy(x, y);
	Copy(M, A);
	
	std::cout << Dot(x,y) << std::endl;
	std::cout << ConjugateDot(x,y) << std::endl;
	
	Scale(x, vector_t::value_type(2));
	Scale(A, matrix_t::value_type(3));
	
	std::cout << "Adding A to M..." << std::endl;
	std::cout << "A:" << std::endl;
	std::cout << A << std::endl;
	std::cout << "M before:" << std::endl;
	std::cout << M << std::endl;
	Add(SubMatrix(A, 0,0,2,2), SubMatrix(M, 1,1,2,2));
	std::cout << "M after:" << std::endl;
	std::cout << M << std::endl;
	
	Add(x, y, vector_t::value_type(2));
	Add(A, M, matrix_t::value_type(2));
	/*
	Rank1Update(A, x, y);
	Rank1Update(A, x);
	Rank2Update(A, x, y);*/
	Mult(A, M, C);
	/*
	std::cout << LargestElementIndex(x) << std::endl;

	A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
	A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;
	A(2,0) = 3; A(2,1) = 8; A(2,2) = 9;

	y[0] = 1; y[1] = 1; y[2] = 1;
	SolveDestructive(A, y);
	std::cout << A << std::endl;
	std::cout << x << std::endl;
	std::cout << y << std::endl;
	*/
	return 0;
}
