#ifndef _RATIONAL_RADICAL1_H_
#define _RATIONAL_RADICAL1_H_

#include "rational.h"
#include <iostream>
#include <cmath>
#include <float.h>
#include <cassert>

/* A number from the field of rationals adjoined with sqrt(n)
 * Non templated
 */

class rational_radical1{
public:
	int n;
	// the value of the number is q + r*sqrt(n)
	rational q, r;
	
	rational_radical1():n(0xDEADBEEF),q(0),r(0){}
	rational_radical1(int val):n(0xDEADBEEF),q(val),r(0){}
	rational_radical1(const rational &val):n(0xDEADBEEF),q(val),r(0){}
	rational_radical1(int root, int val):n(root),q(val),r(0){}
	rational_radical1(int root, const rational &val):n(root),q(val),r(0){}
	rational_radical1(int root, const rational &nonrad, const rational &rad):n(root),q(nonrad),r(rad){}
	rational_radical1(int root, const rational_radical1 &val):n(root),q(val.q),r(val.r){}
	
	rational_radical1& operator=(const rational_radical1 &val){
		n = val.n;
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
	
	int radicand() const{ return n; }
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
	
	rational_radical1& operator+=(const rational_radical1 &val){
		q += val.q;
		if(0 == r.sign()){
			r = val.r;
			n = val.n;
		}else if(0 == val.r.sign()){
			// nothing
		}else{
			assert(radicand() == val.radicand());
			r += val.r;
		}
		return *this;
	}
	rational_radical1& operator-=(const rational_radical1 &val){
		if(0 != val.sign()){ // early exit since sign computation and comparisons use this
			q -= val.q;
			if(0 == r.sign()){
				r -= val.r;
				n = val.n;
			}else if(0 == val.r.sign()){
				// nothing
			}else{
				assert(radicand() == val.radicand());
				r -= val.r;
			}
		}
		return *this;
	}
	rational_radical1& operator*=(const rational_radical1 &val){
		if(0 == r.sign()){
			n = val.n;
			r = q*val.r; q *= val.q;
			return *this;
		}else if(0 == val.r.sign()){
			// nothing
			q *= val.q; r *= val.q;
			return *this;
		}else{
			assert(radicand() == val.radicand());
			// [q1 + r1*sqrt(n)] * [q2 + r2*sqrt(n)]
			// (q1*q2 + r1*r2*n) + (q1*r2 + q2*r1)*sqrt(n)
			return (*this) = rational_radical1(n, q*val.q + r*val.r*n, q*val.r + r*val.q);
		}
	}
	rational_radical1& operator/=(const rational_radical1 &val){
		if(0 == r.sign()){
			n = val.n;
			rational denom(val.q*val.q - val.r*val.r*n);
			return (*this) = rational_radical1(n, (q*val.q)/denom, (- q*val.r)/denom);
		}else if(0 == val.r.sign()){
			// nothing
			q /= val.q; r /= val.q;
			return *this;
		}else{
			assert(radicand() == val.radicand());
			// [q1 + r1*sqrt(n)] / [q2 + r2*sqrt(n)]
			// [(q1*q2 - r1*r2*n) + (r1*q2 - q1*r2)*sqrt(n)] / (q2*q2 - r2*r2*n)
			rational denom(val.q*val.q - val.r*val.r*n);
			return (*this) = rational_radical1(n, (q*val.q - r*val.r*n)/denom, (r*val.q - q*val.r)/denom);
		}

		assert(radicand() == val.radicand());
		rational denom(val.q*val.q - val.r*val.r*n);
		return (*this) = rational_radical1(n, (q*val.q - r*val.r*n)/denom, (r*val.q - q*val.r)/denom);
	}
	
	rational_radical1& operator+=(const rational &val){ q += val; return *this; }
	rational_radical1& operator-=(const rational &val){ q -= val; return *this; }
	rational_radical1& operator*=(const rational &val){ q *= val; r *= val; return *this; }
	rational_radical1& operator/=(const rational &val){ q /= val; r /= val; return *this; }
	
	rational_radical1& operator+=(int val){ q += val; return *this; }
	rational_radical1& operator-=(int val){ q -= val; return *this; }
	rational_radical1& operator*=(int val){ q *= val; r *= val; return *this; }
	rational_radical1& operator/=(int val){ q /= val; r /= val; return *this; }
	
	bool operator==(const rational_radical1& val) const{
		return ((n == val.n) && (q == val.q) && (r == val.r));
	}
	bool operator!=(const rational_radical1& val) const{
		return ((n != val.n) || (q != val.q) || (r != val.r));
	}
	bool operator< (const rational_radical1& val) const;
	bool operator<=(const rational_radical1& val) const;
	bool operator> (const rational_radical1& val) const;
	bool operator>=(const rational_radical1& val) const;
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
	float float_value() const;
	double double_value() const;
	//operator float(){ return float_value(); }
	//operator double(){ return double_value(); }
};

rational_radical1 operator+(const rational_radical1& a, const rational_radical1 &b);
rational_radical1 operator-(const rational_radical1& a, const rational_radical1 &b);
rational_radical1 operator*(const rational_radical1& a, const rational_radical1 &b);
rational_radical1 operator/(const rational_radical1& a, const rational_radical1 &b);
rational_radical1 operator-(const rational_radical1& q);

rational_radical1 operator+(const rational& a, const rational_radical1 &b);
rational_radical1 operator-(const rational& a, const rational_radical1 &b);
rational_radical1 operator*(const rational& a, const rational_radical1 &b);
rational_radical1 operator/(const rational& a, const rational_radical1 &b);

rational_radical1 operator+(int a, const rational_radical1 &b);
rational_radical1 operator-(int a, const rational_radical1 &b);
rational_radical1 operator*(int a, const rational_radical1 &b);
rational_radical1 operator/(int a, const rational_radical1 &b);

rational_radical1 abs(const rational_radical1& val);
int floor(const rational_radical1& q);
int ceil(const rational_radical1& q);
int round(const rational_radical1& q);

std::ostream& operator<<(std::ostream& os, const rational_radical1& q);

#ifdef USING_NUMERIC_TYPE_TRAITS

#include <limits>

template <>
class ScalarTraits<rational_radical1>{
public:
	typedef rational_radical1 value_type;
	
	template <class NumericType>
	static NumericType numeric_value(const value_type &v){ return NumericType(v.double_value()); }
	template <>
	static float numeric_value<float>(const value_type &v){ return v.float_value(); }
	
	template <unsigned int n>
	static const value_type epsilon(){ return rational_radical1(n,ScalarTraits<rational>::epsilon()); }
	template <unsigned int n>
	static const value_type max_value(){ return rational_radical1(n,ScalarTraits<rational>::max_value(), ScalarTraits<rational>::max_value()); }
	template <unsigned int n>
	static const value_type min_value(){ return rational_radical1(n,ScalarTraits<rational>::min_value()); }
};

template <>
class FieldTraits<rational_radical1>{
public:
	typedef rational_radical1 value_type;
	
	static value_type Solve(const value_type &A, const value_type &b, value_type &x){
		x = b/A;
	}
	
	template <unsigned int n>
	static const value_type zero(){ return rational_radical1(n,0); }
	template <unsigned int n>
	static const value_type one(){ return rational_radical1(n,1); }
};

#endif // USING_NUMERIC_TYPE_TRAITS

#endif // _RATIONAL_RADICAL1_H_
