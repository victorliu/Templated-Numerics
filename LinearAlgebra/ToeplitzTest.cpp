#include "ToeplitzMatrix.h"
#include <complex>
#include <iostream>

//typedef std::complex<double> value_type;
typedef double value_type;
typedef ToeplitzMatrix<value_type> Matrix;

void print_matrix(const Matrix &M){
	std::cout << "{";
	for(int i = 0; i < M.size(); ++i){
		std::cout << "{";
		std::cout << M(i,0);
		for(int j = 1; j < M.size(); ++j){
			std::cout << ", " << M(i,j);
		}
		std::cout << "}," << std::endl;
	}
	std::cout << "}" << std::endl;
}

void print_matrix_prim(const value_type *M, size_t N){
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

int main(){
	const int N = 3;
	Matrix M(N);
	for(int i = 0; i < N; ++i){
		M(i,0) = M(0,i) = 1/value_type(i+1);
	}
	M(2,0) = 0.125;
	value_type *inv = new value_type[N*N];
	M.GetInverse(inv, N);
	
	std::cout << "Original matrix:" << std::endl;
	print_matrix(M);
	std::cout << "Inverse matrix:" << std::endl;
	print_matrix_prim(inv, N);
	
	delete [] inv;
	return 0;
}

