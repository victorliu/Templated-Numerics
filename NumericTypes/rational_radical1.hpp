#ifndef _RATIONAL_RADICAL1_HPP_
#define _RATIONAL_RADICAL1_HPP_

#include "rational.h"
#include <iostream>
#include <cmath>
#include <float.h>

/* A number from the field of rationals adjoined with sqrt(n)
 */

template <int n>
class rational_radical1{
public:
	// the value of the number is q + r*sqrt(n)
	rational q, r;
	
	rational_radical1():q(0),r(0){}
	rational_radical1(int val):q(val),r(0){}
	rational_radical1(const rational &val):q(val),r(0){}
	rational_radical1(const rational &nonrad, const rational &rad):q(nonrad),r(rad){}
	rational_radical1(const rational_radical1 &val):q(val.q),r(val.r){}
	
	rational_radical1& operator=(const rational_radical1 &val){
		q = val.q;
		r = val.r;
		return *this;
	}
	rational_radical1& operator=(const rational &val){
		q = val;
		r = 0;
		return *this;
	}
	rational_radical1& operator=(int val){
		q = val;
		r = 0;
		return *this;
	}
	
	const rational& rational_part() const{ return q; }
	const rational& radical_coeff() const{ return r; }
	int sign() const{
		int qs = q.sign();
		int rs = r.sign();
		if(qs == rs){ return qs; }
		if(0 == rs){
			return qs;
		}else if(0 == qs){
			return rs;
		}else{
			rational qr2 = q/r;
			qr2 *= qr2;
			qr2 -= n;
			if(rs > 0){
				return qr2.sign();
			}else{
				return -qr2.sign();
			}
		}
	}
	bool is_nan() const{ return q.is_nan() || r.is_nan(); }
	
	rational_radical1<n>& operator+=(const rational_radical1<n> &val){ q += val.q; r += val.r; return *this; }
	rational_radical1<n>& operator-=(const rational_radical1<n> &val){ q -= val.q; r -= val.r; return *this; }
	rational_radical1<n>& operator*=(const rational_radical1<n> &val){
		// [q1 + r1*sqrt(n)] * [q2 + r2*sqrt(n)]
		// (q1*q2 + r1*r2*n) + (q1*r2 + q2*r1)*sqrt(n)
		return (*this) = rational_radical1<n>(q*val.q + r*val.r*n, q*val.r + r*val.q);
	}
	rational_radical1<n>& operator/=(const rational_radical1<n> &val){
		// [q1 + r1*sqrt(n)] / [q2 + r2*sqrt(n)]
		// [(q1*q2 - r1*r2*n) + (r1*q2 - q1*r2)*sqrt(n)] / (q2*q2 - r2*r2*n)
		rational denom(val.q*val.q - val.r*val.r*n);
		return (*this) = rational_radical1<n>((q*val.q - r*val.r*n)/denom, (r*val.q - q*val.r)/denom);
	}
	/*
	rational_radical1& operator+=(const rational &val);
	rational_radical1& operator-=(const rational &val);
	rational_radical1& operator*=(const rational &val);
	rational_radical1& operator/=(const rational &val);
	
	rational_radical1& operator+=(int val);
	rational_radical1& operator-=(int val);
	rational_radical1& operator*=(int val);
	rational_radical1& operator/=(int val);
	*/
	bool operator==(const rational_radical1& val) const{
		return ((q == val.q) && (r == val.r));
	}
	bool operator!=(const rational_radical1& val) const{
		return ((q != val.q) || (r != val.r));
	}
	bool operator< (const rational_radical1& val) const{
		return ((*this) - val).sign() < 0;
	}
	bool operator<=(const rational_radical1& val) const{
		return ((*this) - val).sign() <= 0;
	}
	bool operator> (const rational_radical1& val) const{
		return ((*this) - val).sign() > 0;
	}
	bool operator>=(const rational_radical1& val) const{
		return ((*this) - val).sign() >= 0;
	}
	/*
	bool operator==(const rational& val) const;
	bool operator!=(const rational& val) const;
	bool operator< (const rational& val) const;
	bool operator<=(const rational& val) const;
	bool operator> (const rational& val) const;
	bool operator>=(const rational& val) const;
	
	bool operator==(int val) const;
	bool operator!=(int val) const;
	bool operator< (int val) const;
	bool operator<=(int val) const;
	bool operator> (int val) const;
	bool operator>=(int val) const;
	*/
	float float_value() const{
		return q.float_value() + r.float_value()*sqrt(n);
	}
	double double_value() const{
		return q.double_value() + r.double_value()*sqrt(n);
	}
	//operator float(){ return float_value(); }
	//operator double(){ return double_value(); }
};

template <int n>
rational_radical1<n> operator+(const rational_radical1<n>& a, const rational_radical1<n> &b){
	rational_radical1<n> ret(a); ret += b; return ret;
}
template <int n>
rational_radical1<n> operator-(const rational_radical1<n>& a, const rational_radical1<n> &b){
	rational_radical1<n> ret(a); ret -= b; return ret;
}
template <int n>
rational_radical1<n> operator*(const rational_radical1<n>& a, const rational_radical1<n> &b){
	rational_radical1<n> ret(a); ret *= b; return ret;
}
template <int n>
rational_radical1<n> operator/(const rational_radical1<n>& a, const rational_radical1<n> &b){
	rational_radical1<n> ret(a); ret /= b; return ret;
}
template <int n>
rational_radical1<n> operator-(const rational_radical1<n>& q){
	return rational_radical1<n>(-q.rational_part(), q.radical_coeff());
}

template <int n>
rational_radical1<n> abs(const rational_radical1<n>& val){
	if(val.sign() >= 0){
		return val;
	}else{
		return -val;
	}
}
template <int n>
int floor(const rational_radical1<n>& q){
	if(q.radical_coeff() == 0){ return floor(q.rational_part()); }
	else{
		// try to use the fact that floor((x+m)/n) == floor((floor(x)+m)/n)
		// so we have
		// n1/d1 + n2/d2*sqrt(n)
		// move to common denominator d, so then we have
		// (n1' + n2'*sqrt(n))/d = (n1' + sqrt(n*n1'*n1'))/d
		// we can use integer square root to get floor of the root
		// This is contingent on the argument to the root being small enough to fit in an int
		int n1 = q.rational_part().num();
		int d1 = q.rational_part().den();
		int n2 = q.radical_coeff().num();
		int d2 = q.radical_coeff().den();
		int sign = q.sign();
		if(sign < 0){ n1 = -n1; n2 = -n2;}
		int g = rational::gcd(d1, d2);
		int dp = d1*(d2/g);
		int n1p = n1*(d2/g);
		int n2p = n2*(d1/g);
		int sqrt_arg = n*n2p*n2p;
		int x = rational::isqrt(sqrt_arg);
		int result = floor(rational(x+n1p, dp));
		return (sign < 0) ? -result : result;
		/*
		double v = q.double_value();
		if(v >= 0){
			return (int)v;
		}else{
			return -(int)(-v);
		}
		*/
	}
}
template <int n>
int ceil(const rational_radical1<n>& q){
	return -floor(-q);
}
template <int n>
int round(const rational_radical1<n>& q){
	if(q.radical_coeff() == 0){ return round(q.rational_part()); }
	else{
		double v = q.double_value();
		if(v >= 0){
			return (int)(v + 0.5);
		}else{
			return -(int)(-v + 0.5);
		}
	}
}

// This will return a negative number if the square root cannot be represented
template <int n>
rational_radical1<n> sqrt(const rational_radical1<n>& val){
	if(val.sign() < 0){
		return val;
	}
	if(0 != val.radical_coeff()){
		return -val;
	}
	rational qn = val.rational_part() / n;
	if(qn * n == val.rational_part()){
		// q must be square rootable
		int n = rational::isqrt(qn.num());
		if(n*n != qn.num()){ return -val; }
		int d = rational::isqrt(qn.den());
		if(d*d != qn.den()){ return -val; }
		return rational_radical<n>(rational(n,d));
	}else{
		// val.rational_part() must be square rootable
		int n = rational::isqrt(q.num());
		if(n*n != q.num()){ return -val; }
		int d = rational::isqrt(q.den());
		if(d*d != q.den()){ return -val; }
		return rational_radical<n>(rational(n,d));
	}
}

template <int n>
std::ostream& operator<<(std::ostream& os, const rational_radical1<n>& q){
	if(q.is_nan()){
		os << "(nan)";
	}else if(q.radical_coeff() == 0){
		os << q.rational_part();
	}else{
		int nu = q.radical_coeff().num();
		int de = q.radical_coeff().den();
		if(q.rational_part() != 0){
			if(nu < 0){
				os << q.rational_part() << '-';
				nu = -nu;
			}else{
				os << q.rational_part() << '+';
			}
		}
		if(1 != nu){
			if(-1 == nu){
				os << '-';
			}else{
				os << nu << '*';
			}
		}
		os << "Sqrt[" << n << ']';
		if(1 != de){
			os << '/' << de;
		}
	}
	return os;
}

#ifdef USING_NUMERIC_TYPE_TRAITS

#include <limits>

template <int n>
class ScalarTraits<rational_radical1<n> >{
public:
	typedef rational_radical1<n> value_type;
	
	template <>
	static double numeric_value<double>(const value_type &v){ return v.double_value(); }
	template <>
	static float numeric_value<float>(const value_type &v){ return v.float_value(); }
	
	static const value_type EPSILON = value_type(ScalarTraits<rational>::EPSILON);
	static const value_type MAX_VALUE = value_type(ScalarTraits<rational>::MAX_VALUE, ScalarTraits<rational>::MAX_VALUE);
	static const value_type MIN_VALUE = value_type(ScalarTraits<rational>::EPSILON);
};

template <int n>
class FieldTraits<rational_radical1<n> >{
public:
	typedef rational_radical1<n> value_type;
	
	static value_type Solve(const value_type &A, const value_type &b, value_type &x){
		x = b/A;
	}
	
	// Must have these:
	static const value_type ZERO = value_type(0);
	static const value_type ONE = value_type(1);
};

#endif // USING_NUMERIC_TYPE_TRAITS

#endif // _RATIONAL_RADICAL1_HPP_
