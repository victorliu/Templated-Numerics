#ifndef _RATIONAL_H_
#define _RATIONAL_H_

#include <iostream>

class rational{
	int n; // carries the sign of the number
	int d; // should always be positive
public:
	rational():n(0),d(1){}
	rational(int numerator):n(numerator),d(1){}
	rational(int numerator, int denominator);
	explicit rational(float value);
	explicit rational(double value);
	rational(const rational &q):n(q.n),d(q.d){}
	
	rational& operator=(const rational &q){ n = q.n; d = q.d; return *this; }
	rational& operator=(int numerator){ n = numerator; d = 1; return *this; }
	
	int num() const{ return n; }
	int den() const{ return d; }
	int sign() const{ return (n==0)?0 : ((n>0)?1:-1); }
	bool is_nan() const{ return (0 == d); }
	
	rational& operator+=(const rational &q);
	rational& operator-=(const rational &q);
	rational& operator*=(const rational &q);
	rational& operator/=(const rational &q);
	
	rational& operator+=(int numerator);
	rational& operator-=(int numerator);
	rational& operator*=(int numerator);
	rational& operator/=(int numerator);
	
	const rational& operator++();
	const rational& operator--();
	
	bool operator==(const rational& q) const;
	bool operator!=(const rational& q) const;
	bool operator< (const rational& q) const;
	bool operator<=(const rational& q) const;
	bool operator> (const rational& q) const;
	bool operator>=(const rational& q) const;
	
	bool operator==(int numerator) const;
	bool operator!=(int numerator) const;
	bool operator< (int numerator) const;
	bool operator<=(int numerator) const;
	bool operator> (int numerator) const;
	bool operator>=(int numerator) const;
	
	float float_value() const;
	double double_value() const;
	//operator float(){ return float_value(); }
	//operator double(){ return double_value(); }
public:
	static unsigned int gcd(unsigned int a, unsigned int b);
	static unsigned int lcm(unsigned int a, unsigned int b);
	static int isqrt(int a);
private:
	static double double_to_rational(double x, int &numerator, int &denominator);
};

rational operator+(const rational& a, const rational &b);
rational operator-(const rational& a, const rational &b);
rational operator*(const rational& a, const rational &b);
rational operator/(const rational& a, const rational &b);
rational operator-(const rational& q);

rational abs(const rational& q);
int floor(const rational& q);
int ceil(const rational& q);
int round(const rational& q);

std::istream& operator>>(std::istream& is, rational& q);
std::ostream& operator<<(std::ostream& os, const rational& q);

#ifdef USING_NUMERIC_TYPE_TRAITS

#include <limits>

template <>
class ScalarTraits<rational>{
public:
	typedef rational value_type;
	
	template <class NumericType>
	static NumericType numeric_value(const value_type &v){ return NumericType(v.double_value()); }
	
	static const value_type epsilon(){ return rational(1,std::numeric_limits<int>::max()); }
	static const value_type max_value(){ return rational(std::numeric_limits<int>::max()); }
	static const value_type min_value(){ return rational(1,std::numeric_limits<int>::max()); }
};

template <>
float ScalarTraits<rational>::numeric_value<float>(const value_type &v);

template <>
class FieldTraits<rational>{
public:
	typedef rational value_type;
	
	static void Solve(const value_type &A, const value_type &b, value_type &x){
		x = b/A;
	}
	
	// Must have these:
	static const value_type zero(){ return rational(0); }
	static const value_type one(){ return rational(1); }
};

#endif // USING_NUMERIC_TYPE_TRAITS

#endif // _RATIONAL_H_
