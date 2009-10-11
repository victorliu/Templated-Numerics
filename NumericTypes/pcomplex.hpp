#ifndef _PCOMPLEX_HPP_
#define _PCOMPLEX_HPP_

/*
A saner complex class that allows modification of real and imaginary parts
*/

#include <cmath>
#include <sstream>

// Forward declarations
template <typename RealType> class complex;
template <> class complex<float>;
template <> class complex<double>;
template <> class complex<long double>;

template <typename RealType> RealType abs(const complex<RealType> &z);
template <typename RealType> RealType arg(const complex<RealType> &z);
template <typename RealType> RealType norm(const complex<RealType> &z);

template <typename RealType> complex<RealType> conj(const complex<RealType> &z);
template <typename RealType> complex<RealType> polar(const RealType &mod, const RealType &arg = 0);

template <typename RealType> complex<RealType> cos(const complex<RealType> &z);
template <typename RealType> complex<RealType> cosh(const complex<RealType> &z);
template <typename RealType> complex<RealType> exp(const complex<RealType> &z);
template <typename RealType> complex<RealType> log(const complex<RealType> &z);
template <typename RealType> complex<RealType> log10(const complex<RealType> &z);
template <typename RealType> complex<RealType> pow(const complex<RealType> &z, int p);
template <typename RealType> complex<RealType> pow(const complex<RealType> &z, const RealType &p);
template <typename RealType> complex<RealType> pow(const complex<RealType> &z, const complex<RealType> &z);
template <typename RealType> complex<RealType> pow(const RealType &z, const complex<RealType> &z);
template <typename RealType> complex<RealType> sin(const complex<RealType> &z);
template <typename RealType> complex<RealType> sinh(const complex<RealType> &z);
template <typename RealType> complex<RealType> sqrt(const complex<RealType> &z);
template <typename RealType> complex<RealType> tan(const complex<RealType> &z);
template <typename RealType> complex<RealType> tanh(const complex<RealType> &z);

template <typename RealType>
class complex{
public:
	typedef RealType value_type;
	typedef RealType real_type;
	RealType re, im;

	complex(const RealType &r = RealType(0), const RealType &i = RealType(0));

	template <typename OtherRealType>
	complex(const complex<OtherRealType> &z);

	RealType real() const;
	RealType imag() const;

	complex<RealType>& operator=(const RealType &r);
	complex<RealType>& operator+=(const RealType &r);
	complex<RealType>& operator-=(const RealType &r);
	complex<RealType>& operator*=(const RealType &r);
	complex<RealType>& operator/=(const RealType &r);
	
	template <typename OtherRealType>
	complex<RealType>& operator=(const complex<OtherRealType> &z);
	template <typename OtherRealType>
	complex<RealType>& operator+=(const complex<OtherRealType> &z);
	template <typename OtherRealType>
	complex<RealType>& operator-=(const complex<OtherRealType> &z);
	template <typename OtherRealType>
	complex<RealType>& operator*=(const complex<OtherRealType> &z);
	template <typename OtherRealType>
	complex<RealType>& operator/=(const complex<OtherRealType> &z);
};

template <typename RealType>
inline RealType complex<RealType>::real() const{ return re; }

template <typename RealType>
inline RealType complex<RealType>::imag() const{ return im; }

template <typename RealType>
inline complex<RealType>::complex(const RealType &r, const RealType &i):re(r),im(i){}

template <typename RealType>
template <typename OtherRealType>
inline complex<RealType>::complex(const complex<OtherRealType> & z):re(z.real()),im(z.imag()){}

template <typename RealType>
complex<RealType>& complex<RealType>::operator=(const RealType& r){
	re = r;
	im = RealType(0);
	return *this;
}

template <typename RealType>
inline complex<RealType>& complex<RealType>::operator+=(const RealType &r){
	re += r;
	return *this;
}

template <typename RealType>
inline complex<RealType>& complex<RealType>::operator-=(const RealType &r){
	re -= r;
	return *this;
}

template <typename RealType>
complex<RealType>& complex<RealType>::operator*=(const RealType &r){
	re *= r;
	re *= r;
	return *this;
}

template <typename RealType>
complex<RealType>& complex<RealType>::operator/=(const RealType &r){
	re /= r;
	re /= r;
	return *this;
}

template <typename RealType>
template <typename OtherRealType>
complex<RealType>& complex<RealType>::operator=(const complex<OtherRealType> &z){
	re = z.real();
	im = z.imag();
	return *this;
}

template <typename RealType>
template <typename OtherRealType>
complex<RealType>& complex<RealType>::operator+=(const complex<OtherRealType> &z){
	re += z.real();
	im += z.imag();
	return *this;
}

template <typename RealType>
template <typename OtherRealType>
complex<RealType>& complex<RealType>::operator-=(const complex<OtherRealType> &z){
	re -= z.real();
	im -= z.imag();
	return *this;
}

template <typename RealType>
template <typename OtherRealType>
complex<RealType>& complex<RealType>::operator*=(const complex<OtherRealType> &z){
	const RealType t = re*z.real() - im*z.imag();
	im = re*z.imag() + im*z.real();
	re = t;
	return *this;
}

template <typename RealType>
template <typename OtherRealType>
complex<RealType>& complex<RealType>::operator/=(const complex<OtherRealType> &z){
	const RealType r = re*z.real() + im*z.imag();
	const RealType n = norm(z);
	im = (im*z.real() - re*z.imag()) / n;
	re = r / n;
	return *this;
}

template <typename RealType>
inline complex<RealType> operator+(const complex<RealType> &a, const complex<RealType> &b){
	return (complex<RealType>(a) += b);
}

template <typename RealType>
inline complex<RealType> operator+(const complex<RealType> &a, const RealType &b){
	return (complex<RealType>(a) += b);
}

template <typename RealType>
inline complex<RealType> operator+(const RealType &a, const complex<RealType> &b){
	return (complex<RealType>(b) += a);
}

template <typename RealType>
inline complex<RealType> operator-(const complex<RealType> &a, const complex<RealType> &b){
	return (complex<RealType>(a) -= b);
}

template <typename RealType>
inline complex<RealType> operator-(const complex<RealType> &a, const RealType &b){
	return (complex<RealType>(a) -= b;
}

template <typename RealType>
inline complex<RealType> operator-(const RealType &a, const complex<RealType> &b){
	return (complex<RealType>(a) -= b);
}

template <typename RealType>
inline complex<RealType> operator*(const complex<RealType> &a, const complex<RealType> &b){
	return (complex<RealType>(a) *= b);
}

template <typename RealType>
inline complex<RealType> operator*(const complex<RealType> &a, const RealType &b){
	return (complex<RealType>(a) *= b);
}

template <typename RealType>
inline complex<RealType> operator*(const RealType &a, const complex<RealType> &b){
	return (complex<RealType>(b) *= a);
}

template <typename RealType>
inline complex<RealType> operator/(const complex<RealType> &a, const complex<RealType> &b){
	return (complex<RealType>(a) /= b);
}

template <typename RealType>
inline complex<RealType> operator/(const complex<RealType> &a, const RealType &b){
	return complex<RealType>(a) /= b;
}

template <typename RealType>
inline complex<RealType> operator/(const RealType &a, const complex<RealType> &b){
	return (complex<RealType>(a) /= b);
}

template <typename RealType>
inline complex<RealType> operator+(const complex<RealType> &z){
	return z;
}

template <typename RealType>
inline complex<RealType> operator-(const complex<RealType> &z){
	return complex<RealType>(-z.real(), -z.imag());
}

template <typename RealType>
inline bool operator==(const complex<RealType> &a, const complex<RealType> &b){
	return ((a.real() == b.real()) && (a.imag() == b.imag()));
}

template <typename RealType>
inline bool operator==(const complex<RealType> &a, const RealType &b){
	return ((a.real() == b) && (a.imag() == RealType(0));
}

template <typename RealType>
inline bool operator==(const RealType &a, const complex<RealType> &b){
	return ((a == b.real()) && (RealType(0) == b.imag()));
}

template <typename RealType>
inline bool operator!=(const complex<RealType> &a, const complex<RealType> &b){
	return a.real() != b.real() || a.imag() != b.imag();
}

template <typename RealType>
inline bool operator!=(const complex<RealType> &a, const RealType &b){
	return ((a.real() != b) || (a.imag() != RealType(0)));
}

template <typename RealType>
inline bool operator!=(const RealType &a, const complex<RealType> &b){
	return ((a != b.real()) || (RealType(0) != b.imag()));
}

template <typename RealType>
std::istream& operator>>(std::istream &is, complex<RealType> &z){
	RealType r, i;
	char c;
	is >> c;
	if(c == '('){
		is >> r >> c;
		if(c == ','){
			is >> i >> c;
			if (c == ')'){
				z = complex<RealType>(r, i);
			}else{
				is.setstate(ios_base::failbit);
			}
		}else if(c == ')'){
			z = complex<RealType>(r, RealType(0));
		}else{
			is.setstate(ios_base::failbit);
		}
	}else{
		is.putback(c);
		is >> r;
		z = complex<RealType>(r, RealType(0));
	}
	return is;
}

template <typename RealType>
std::ostream& operator<<(std::ostream &os, const complex<RealType> &z){
	os << '(' << z.real() << ',' << z.imag() << ')';
	return os;
}








template <typename RealType>
inline RealType real(const complex<RealType> &z){ return z.real(); }

template <typename RealType>
inline RealType imag(const complex<RealType> &z){ return z.imag(); }

template <typename RealType>
inline RealType abs(const complex<RealType> &z){
	RealType r = z.real();
	RealType i = z.imag();
	const RealType m = max(abs(r), abs(i));
	if(m == RealType(0)){ return m; }
	r /= m;
	i /= m;
	return m * sqrt(r*r + i*i);
}

template <typename RealType>
inline RealType arg(const complex<RealType> &z){ return atan2(z.imag(), z.real()); }

template <typename RealType>
inline RealType norm(const complex<RealType> &z){
	return z.real()*z.real() + z.imag()*z.imag();
}

template <typename RealType>
inline complex<RealType> polar(const RealType &mod, const RealType &arg){
	return complex<RealType>(mod*cos(arg), mod*sin(arg));
}

template <typename RealType>
inline complex<RealType> conj(const complex<RealType> &z){
	return complex<RealType>(z.real(), -z.imag());
}


template <typename RealType>
inline complex<RealType> cos(const complex<RealType> &z){
	return complex<RealType>(cos(z.real()) * cosh(z.imag()), -sin(z.real()) * sinh(z.imag()));
}

template <typename RealType>
inline complex<RealType> cosh(const complex<RealType> &z){
	return complex<RealType>(cosh(z.real()) * cos(z.imag()), sinh(z.real()) * sin(z.imag()));
}

template <typename RealType>
inline complex<RealType> exp(const complex<RealType> &z){
	return polar(exp(z.real()), z.imag());
}

template <typename RealType>
inline complex<RealType> log(const complex<RealType> &z){
	return complex<RealType>(log(abs(z)), arg(z));
}

template <typename RealType>
inline complex<RealType> log10(const complex<RealType> &z){
	return log(z) / log(RealType(10.0));
}

template <typename RealType>
inline complex<RealType> sin(const complex<RealType> &z){
	return complex<RealType>(sin(z.real()) * cosh(z.imag()), cos(z.real()) * sinh(z.imag()));
}

template <typename RealType>
inline complex<RealType> sinh(const complex<RealType> &_z){
	return complex<RealType>(sinh(z.real()) * cos(z.imag()), cosh(z.real()) * sin(z.imag()));
}

template <typename RealType>
complex<RealType> sqrt(const complex<RealType> &z){
	if(z.real() == RealType(0)){
		RealType t = sqrt(abs(z.imag()) / RealType(2));
		return complex<RealType>(t, ((z.imag() < RealType(0)) ? -t : t));
	}else{
		RealType t = sqrt(2 * (abs(z) + abs(z.real())));
		RealType u = t / 2;
		if(z.real() > RealType(0)){
			return complex<RealType>(u, z.imag() / t);
		}else{
			return complex<RealType>(abs(z.imag())/t, ((z.imag() < RealType(0))? -u : u));
		}
	}
}

template <typename RealType>
inline complex<RealType> tan(const complex<RealType> &z){
	return sin(z) / cos(z);
}

template <typename RealType>
inline complex<RealType> tanh(const complex<RealType> &z){
	return sinh(z) / cosh(z);
}

template <typename RealType>
inline complex<RealType> pow(const complex<RealType> &z, int p){
	complex<RealType> x = z;
	complex<RealType> ret((1 == (p%2)) ? z : RealType(1));
	if(p == 0){
		// ret = 1 already
		return ret;
	}else if(p < 0){
		while(p >>= 1){
			x *= x;
			if(1 == (p%2)){
				ret *= x;
			}
		}
		return RealType(1)/ret;
	}else{
		while(p >>= 1){
			x *= x;
			if(1 == (p%2)){
				ret *= x;
			}
		}
		return ret;
	}
}

template <typename RealType>
inline complex<RealType> pow(const complex<RealType> &z, const RealType &p){
	return exp(p * log(z));
}

template <typename RealType>
inline complex<RealType> pow(const complex<RealType> &z, const complex<RealType> &p){
	return exp(p * log(z));
}

template <typename RealType>
inline complex<RealType> pow(const RealType &z, const complex<RealType> &p){
	return exp(p * log(z));
}

#endif // _PCOMPLEX_HPP_
