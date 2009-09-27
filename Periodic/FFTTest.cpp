#include "FFT.hpp"
#include <complex>
#include <iostream>

typedef std::complex<double> complex_t;

int main(){
	const size_t N = 14;
	complex_t *data = new complex_t[N];
	
	for(size_t i = 0; i < N; ++i){
		data[i] = i;
	}
	
	std::cout << "Original:" << std::endl;
	for(size_t i = 0; i < N; ++i){
		std::cout << data[i] << std::endl;
	}
	FFT::Transform(N, data, data, FFT::FORWARD);
	std::cout << "Transformed:" << std::endl;
	for(size_t i = 0; i < N; ++i){
		std::cout << data[i] << std::endl;
	}
	FFT::Transform(N, data, data, FFT::INVERSE);
	std::cout << "Inversed:" << std::endl;
	for(size_t i = 0; i < N; ++i){
		std::cout << data[i] << std::endl;
	}
	
	delete [] data;
	return 0;
}
