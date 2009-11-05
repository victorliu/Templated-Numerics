#include "rational_radical1.h"
#include <iostream>
#include <cmath>
#include <float.h>

rational_radical1 operator+(const rational_radical1& a, const rational_radical1 &b){
	rational_radical1 ret(a); ret += b; return ret;
}
rational_radical1 operator-(const rational_radical1& a, const rational_radical1 &b){
	rational_radical1 ret(a); ret -= b; return ret;
}
rational_radical1 operator*(const rational_radical1& a, const rational_radical1 &b){
	rational_radical1 ret(a); ret *= b; return ret;
}
rational_radical1 operator/(const rational_radical1& a, const rational_radical1 &b){
	rational_radical1 ret(a); ret /= b; return ret;
}
rational_radical1 operator-(const rational_radical1& q){
	return rational_radical1(q.radicand(), -q.rational_part(), q.radical_coeff());
}

rational_radical1 operator+(const rational& a, const rational_radical1 &b){
	return (rational_radical1(b) += a);
}
rational_radical1 operator-(const rational& a, const rational_radical1 &b){
	return (rational_radical1(b.radicand(), a) -= b);
}
rational_radical1 operator*(const rational& a, const rational_radical1 &b){
	return (rational_radical1(b) *= a);
}
rational_radical1 operator/(const rational& a, const rational_radical1 &b){
	return (rational_radical1(b.radicand(), a) /= b);
}

rational_radical1 operator+(int a, const rational_radical1 &b){
	return (rational_radical1(b) += a);
}
rational_radical1 operator-(int a, const rational_radical1 &b){
	return (rational_radical1(b.radicand(), a) -= b);
}
rational_radical1 operator*(int a, const rational_radical1 &b){
	return (rational_radical1(b) *= a);
}
rational_radical1 operator/(int a, const rational_radical1 &b){
	return (rational_radical1(b.radicand(), a) /= b);
}

float rational_radical1::float_value() const{
	return q.float_value() + r.float_value()*sqrt((float)n);
}
double rational_radical1::double_value() const{
	return q.double_value() + r.double_value()*sqrt((double)n);
}


bool rational_radical1::operator< (const rational_radical1& val) const{
	return ((*this) - val).sign() < 0;
}
bool rational_radical1::operator<=(const rational_radical1& val) const{
	return ((*this) - val).sign() <= 0;
}
bool rational_radical1::operator> (const rational_radical1& val) const{
	return ((*this) - val).sign() > 0;
}
bool rational_radical1::operator>=(const rational_radical1& val) const{
	return ((*this) - val).sign() >= 0;
}

rational_radical1 abs(const rational_radical1& val){
	if(val.sign() >= 0){
		return val;
	}else{
		return -val;
	}
}
int floor(const rational_radical1& q){
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
		int sqrt_arg = q.radicand()*n2p*n2p;
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
int ceil(const rational_radical1& q){
	return -floor(-q);
}
int round(const rational_radical1& q){
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

std::ostream& operator<<(std::ostream& os, const rational_radical1& q){
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
		os << "Sqrt[" << q.radicand() << ']';
		if(1 != de){
			os << '/' << de;
		}
	}
	return os;
}

