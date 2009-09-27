#ifndef _FFT_H_
#define _FFT_H_

#include <cstdlib>
#include <cmath>
#include <vector>

// Most of the code is taken from Kiss FFT, ported by Victor Liu

/*
Copyright (c) 2003-2006 Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

namespace FFT{

enum Direction{
	FORWARD,
	INVERSE
};

template <class ComplexType>
class Plan{
	Direction dir;
	std::vector<size_t> factors;
	std::vector<ComplexType> twiddles;
public:
	Plan(size_t length, Direction Dir = FORWARD):dir(Dir){
		size_t n = length;
		
		// Compute the twiddles
		twiddles.reserve(n);
		for(size_t i = 0; i < n; ++i){
			typename ComplexType::value_type phase = -2*M_PI*i/n;
			if(INVERSE == dir){ phase = -phase; }
			twiddles.push_back(ComplexType(cos(phase), sin(phase)));
		}
		
		// Factorization
		
		factors.push_back(n);
		// Pull out all factors of 4
		while(0 == (n&3)){
			factors.push_back(4);
			n >>= 2;
			factors.push_back(n);
		}
		// Pull out factor of 2
		if(0 == (n&1)){
			factors.push_back(2);
			n >>= 1;
			factors.push_back(n);
		}

		size_t rootn = (size_t)sqrt((double)n+0.5);
		for(size_t i = 3; i <= rootn; i += 2){
			size_t q;
			while((q = (n/i))*i == n){
				factors.push_back(i);
				n = q;
				factors.push_back(n);
				rootn = (size_t)sqrt((double)n+0.5); // update the limit
			}
		}
		if(n > 1){ // should never get here
			factors.push_back(n);
			factors.push_back(1);
		}
		/*
		for(int i = 0; i < factors.size(); ++i){
			std::cout << factors[i] << std::endl;
		}
		*/
	}
	size_t size() const{ return factors[0]; }
	Direction direction() const{ return dir; }
	size_t operator[](size_t i) const{ return factors[i+1]; }
	const ComplexType& operator()(size_t i) const{ return twiddles[i]; }
};


template <class ComplexType>
static void TransformButterfly2(const Plan<ComplexType> &plan, ComplexType *out, const size_t stride, size_t m){
	ComplexType *out2;
	const ComplexType *tw1 = &(plan(0));
	ComplexType t;
	out2 = out + m;
	do{
		out[0] /= 2;
		out2[0] /= 2;

		t = out2[0] * tw1[0];
		tw1 += stride;
		out2[0] = out[0] - t;
		out[0] += t;
		++out2;
		++out;
	}while(--m);
}

template <class ComplexType>
static void TransformButterfly4(const Plan<ComplexType> &plan, ComplexType *out, const size_t stride, size_t m){
	const ComplexType *tw1,*tw2,*tw3;
	ComplexType scratch[6];
	size_t k=m;
	const size_t m2=2*m;
	const size_t m3=3*m;

	tw3 = tw2 = tw1 = &(plan(0));

	do{
		out[0] /= 4;
		out[m] /= 4;
		out[m2] /= 4;
		out[m3] /= 4;

		scratch[0] = out[m] * (*tw1);
		scratch[1] = out[m2] * (*tw2);
		scratch[2] = out[m3] * (*tw3);

		scratch[5] = *out - scratch[1];
		out[0] += scratch[1];
		scratch[3] = scratch[0] + scratch[2];
		scratch[4] = scratch[0] - scratch[2];
		out[m2] = *out - scratch[3];
		tw1 += stride;
		tw2 += stride*2;
		tw3 += stride*3;
		out[0] += scratch[3];

		if(INVERSE == plan.direction()){
			out[m] = ComplexType(scratch[5].real() - scratch[4].imag(), scratch[5].imag() + scratch[4].real());
			out[m3] = ComplexType(scratch[5].real() + scratch[4].imag(), scratch[5].imag() - scratch[4].real());
		}else{
			out[m] = ComplexType(scratch[5].real() + scratch[4].imag(), scratch[5].imag() - scratch[4].real());
			out[m3] = ComplexType(scratch[5].real() - scratch[4].imag(), scratch[5].imag() + scratch[4].real());
		}
		++out;
	}while(--k);
}

template <class ComplexType>
static void TransformButterfly3(const Plan<ComplexType> &plan, ComplexType *out, const size_t stride, size_t m){
	size_t k=m;
	const size_t m2 = 2*m;
	const ComplexType *tw1,*tw2;
	ComplexType scratch[5];
	ComplexType epi3 = plan(stride*m);

	tw1=tw2= &(plan(0));

	do{
		out[0] /= 3;
		out[m] /= 3;
		out[m2] /= 3;

		scratch[1] = out[m] * (*tw1);
		scratch[2] = out[m2] * (*tw2);

		scratch[3] = scratch[1] + scratch[2];
		scratch[0] = scratch[1] - scratch[2];
		tw1 += stride;
		tw2 += stride*2;

		out[m] = ComplexType(out[0].real() - scratch[3].real()/2, out[0].imag() - scratch[3].imag()/2);

		scratch[0] *= epi3.imag();

		out[0] += scratch[3];

		out[m2] = ComplexType(out[m].real() + scratch[0].imag(), out[m].imag() - scratch[0].real());

		out[m] = ComplexType(out[m].real() - scratch[0].imag(), out[m].imag() + scratch[0].real());

		++out;
	}while(--k);
}

template <class ComplexType>
static void TransformButterfly5(const Plan<ComplexType> &plan, ComplexType *out, const size_t stride, size_t m){
	ComplexType *out0,*out1,*out2,*out3,*out4;
	ComplexType scratch[13];
	const ComplexType *twiddles = &(plan(0));
	const ComplexType *tw;
	ComplexType ya,yb;
	ya = twiddles[stride*m];
	yb = twiddles[stride*2*m];

	out0=out;
	out1=out0+m;
	out2=out0+2*m;
	out3=out0+3*m;
	out4=out0+4*m;

	tw= &(plan(0));
	for(size_t u=0; u<m; ++u ) {
		out0[0] /= 5;
		out1[0] /= 5;
		out2[0] /= 5;
		out3[0] /= 5;
		out4[0] /= 5;
		scratch[0] = out0[0];

		scratch[1] = out1[0] * tw[1*u*stride];
		scratch[2] = out2[0] * tw[2*u*stride];
		scratch[3] = out3[0] * tw[3*u*stride];
		scratch[4] = out4[0] * tw[4*u*stride];

		scratch[7] = scratch[1] + scratch[4];
		scratch[10] = scratch[1] - scratch[4];
		scratch[8] = scratch[2] + scratch[3];
		scratch[9] = scratch[2] - scratch[3];

		out0[0] += scratch[7] + scratch[8];

		scratch[5] = ComplexType(scratch[0].real() + (scratch[7].real() * ya.real()) + (scratch[8].real() * yb.real()),
		                         scratch[0].imag() + (scratch[7].imag() * ya.real()) + (scratch[8].imag() * yb.real()));
		
		scratch[6] = ComplexType((scratch[10].imag() * ya.imag()) + (scratch[9].imag() * yb.imag()), -(scratch[10].real() * ya.imag()) - (scratch[9].real() * yb.imag()));

		out1[0] = scratch[5] - scratch[6];
		out4[0] = scratch[5] + scratch[6];

		scratch[11] = ComplexType(scratch[0].real() + (scratch[7].real()*yb.real()) + (scratch[8].real()*ya.real()), scratch[0].imag() + (scratch[7].imag()*yb.real()) + (scratch[8].imag()*ya.real()));
		scratch[12] = ComplexType(- (scratch[10].imag()*yb.imag()) + (scratch[9].imag()*ya.imag()), (scratch[10].real()*yb.imag()) - (scratch[9].real()*ya.imag()));

		out2[0] = scratch[11] + scratch[12];
		out3[0] = scratch[11] - scratch[12];

		++out0;++out1;++out2;++out3;++out4;
	}
}

template <class ComplexType>
static void TransformButterflyGeneric(const Plan<ComplexType> &plan, ComplexType *out, const size_t stride, size_t m, size_t p){
	const ComplexType *twiddles = &(plan(0));
	int Norig = plan.size();

	std::vector<ComplexType> scratchbuf(p);

	for(size_t u=0; u < m; ++u){
		size_t k = u;
		for(size_t q1 = 0; q1 < p; ++q1){
			scratchbuf[q1] = out[k];
			scratchbuf[q1] /= p;
			k += m;
		}

		k = u;
		for(size_t q1 = 0; q1 < p; ++q1){
			size_t twidx = 0;
			out[ k ] = scratchbuf[0];
			for(size_t q = 1; q < p; ++q){
				twidx += stride * k;
				if(twidx >= Norig){ twidx -= Norig; }
				out[k] += scratchbuf[q] * twiddles[twidx];
			}
			k += m;
		}
	}
}

template <class ComplexType>
static void TransformStep(const Plan<ComplexType> &plan, const ComplexType *in, ComplexType *out, size_t stride, size_t which_factor){
	ComplexType *out_beg = out;
	const size_t p = plan[which_factor++]; // the radix
	const size_t m = plan[which_factor++]; // stage's fft length/p
	const ComplexType * out_end = out + p*m;

	if(1 == m){
		do{
			*out = *in;
			in += stride;
		}while(++out != out_end );
	}else{
		do{
			// recursive call:
			// DFT of size m*p performed by doing
			// p instances of smaller DFTs of size m, 
			// each one takes a decimated version of the input
			TransformStep(plan, in, out, stride*p, which_factor);
			in += stride;
		}while((out += m) != out_end);
	}

	out = out_beg;

	// recombine the p smaller DFTs 
	switch(p){
		case 2:  TransformButterfly2      (plan, out, stride, m   ); break;
		case 3:  TransformButterfly3      (plan, out, stride, m   ); break; 
		case 4:  TransformButterfly4      (plan, out, stride, m   ); break;
		case 5:  TransformButterfly5      (plan, out, stride, m   ); break; 
		default: TransformButterflyGeneric(plan, out, stride, m, p); break;
	}
}

template <class ComplexType>
void Transform(const Plan<ComplexType> &plan, const ComplexType *in, ComplexType *out){
	ComplexType *out2 = out;
	if(in == out){
		out2 = new ComplexType[plan.size()];
	}
	
	TransformStep(plan, in, out2, 1, 0);
	
	if(in == out){
		memcpy(out, out2, sizeof(ComplexType)*plan.size());
		delete [] out2;
	}
}
template <class ComplexType>
void Transform(size_t n, const ComplexType *in, ComplexType *out, Direction dir = FORWARD){
	Plan<ComplexType> plan(n, dir);
	Transform(plan, in, out);
}

}; // namesapce FFT

#endif // _FFT_H_
