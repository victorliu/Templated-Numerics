#ifndef _COMPLEX_IO_HPP_
#define _COMPLEX_IO_HPP_

#include <iostream>
#include <complex>
#include <limits>

template <class T>
std::ostream&
	operator<<(
		std::ostream& ostr,
		const std::complex<T>& z
	)
{
	if(0 == z.imag() || std::abs(z.real()*std::numeric_limits<T>::epsilon()) > std::abs(z.imag())){
		ostr << z.real();
	}else{
		if(z.imag() < 0){
			ostr << z.real() << "-I*" << (-z.imag());
		}else{
			ostr << z.real() << "+I*" << z.imag();
		}
	}
	return ostr;
}

#endif // _COMPLEX_IO_HPP_
