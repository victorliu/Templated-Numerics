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
void print_matrix_inverse(const Matrix::InverseMatrix &M){
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
	const int N = 4;
	Matrix M(N);
	for(int i = 0; i < N; ++i){
		M(i,0) = M(0,i) = 1/value_type(i+1);
	}
	M(2,0) = 0.125;
	value_type *invraw = new value_type[N*N];
	value_type *orig = new value_type[N*N]; M.Fill(orig, N);
	Matrix::InverseMatrix inv(M);
	inv.Fill(invraw, N);
	
	std::cout << "Original matrix:" << std::endl;
	print_matrix(M);
	print_matrix_prim(orig,N);
	
	std::cout << "Inverse matrix:" << std::endl;
	print_matrix_prim(invraw, N);
	print_matrix_inverse(inv);
	
	value_type *scratch = new value_type[N*N];memset(scratch,0,N*N*sizeof(value_type));
	inv.MultAdd(N, orig, N, scratch, N);
	std::cout << "Inverse*Original:" << std::endl;
	print_matrix_prim(scratch, N);
	
	delete [] invraw;
	delete [] scratch;
	return 0;
}

