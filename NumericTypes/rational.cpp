#include "rational.h"
#include <cmath>
#include <cstdlib>

rational::rational(int numerator, int denominator):n(numerator),d(denominator){
	if(d < 0){
		n = -n;
		d = -d;
	}
	if(0 == n){ d = 1; return; }
	int g = gcd(std::abs((int)n), d);
	n /= g; d /= g;
}
rational::rational(float value){
	double_to_rational((double)value, n, d);
}
rational::rational(double value){
	double_to_rational(value, n, d);
}

rational& rational::operator+=(const rational &q){
	if(is_nan()){
		if(q.is_nan()){ // addition of infinities
			n += q.n;
			if(n < 0){ n = -1; }
			else if(n > 0){ n = 1; }
		}
		return *this;
	}else if(q.is_nan()){
		n = q.n;
		d = q.d;
		return *this;
	}
	// stolen from boost::rational:
	// We have to compute a/b + c/d, where gcd(a,b)=1 and gcd(c,d)=1.
	// Let g = gcd(b,d), and b1 = b/g, d1 = d/g. Then gcd(b1,d1)=1
	//
	// The result is (a*d1 + c*b1) / (b1*d1*g).
	// Now we have to normalize this ratio.
	// Let's assume h | gcd((a*d1 + c*b1), (b1*d1*g)), and h > 1
	// If h | b1 then gcd(h,d1)=1 and hence h|(a*d1+c*b1) => h|a.
	// But since gcd(a,b1)=1 we have h=1.
	// Similarly h|d1 leads to h=1.
	// So we have that h | gcd((a*d1 + c*b1) , (b1*d1*g)) => h|g
	// Finally we have gcd((a*d1 + c*b1), (b1*d1*g)) = gcd((a*d1 + c*b1), g)
	// Which proves that instead of normalizing the result, it is better to
	// divide num and den by gcd((a*d1 + c*b1), g)
	
	int g = gcd(d, q.d);
	d /= g;  // = b1 from the calculations above
	n = n * (q.d / g) + q.n * d;
	g = gcd(abs(n), abs(g));
	n /= g;
	d *= q.d/g;

	return *this;
}
rational& rational::operator-=(const rational &q){
	if(is_nan()){
		if(q.is_nan()){ // addition of infinities
			n -= q.n;
			if(n < 0){ n = -1; }
			else if(n > 0){ n = 1; }
		}
		return *this;
	}else if(q.is_nan()){
		n = q.n;
		d = q.d;
		return *this;
	}
	int g = gcd(d, q.d);
	d /= g;  // = b1 from the calculations above
	n = n * (q.d / g) - q.n * d;
	g = gcd(abs(n), g);
	n /= g;
	d *= q.d/g;

	return *this;
}
rational& rational::operator*=(const rational &q){
	if(is_nan()){
		if(q.is_nan()){ // multiplication of infinities
			n *= q.sign();
		}else{
			n *= q.sign();
		}
		return *this;
	}else if(q.is_nan()){
		n = q.n * sign();
		d = q.d;
		return *this;
	}
	
	int g1 = gcd(abs(n), q.d);
	int g2 = gcd(abs(q.n), d);
	n = (n/g1) * (q.n/g2);
	d = (d/g2) * (q.d/g1);
	return *this;
}
rational& rational::operator/=(const rational &q){
	if(is_nan()){
		if(q.is_nan()){ // division of infinities
			n = 0;
			d = 0;
		}else{
			n *= q.sign();
		}
		return *this;
	}else if(q.is_nan()){
		n = 0;
		d = 0;
		return *this;
	}
	
	if(q == 0){
		n = sign();
		d = 0;
		return *this;
	}else if(0 == n){
		return *this;
	}
	
	int g1 = gcd(abs(n), abs(q.n));
	int g2 = gcd(q.d, d);
	n = (n/g1) * (q.d/g2);
	d = (d/g2) * (q.n/g1);

	if(d < 0){
		n = -n;
		d = -d;
	}
	return *this;
}

rational& rational::operator+=(int numerator){ return ((*this) += rational(numerator)); }
rational& rational::operator-=(int numerator){ return ((*this) -= rational(numerator)); }
rational& rational::operator*=(int numerator){ return ((*this) *= rational(numerator)); }
rational& rational::operator/=(int numerator){ return ((*this) /= rational(numerator)); }

const rational& rational::operator++(){ n += d; return *this; }
const rational& rational::operator--(){ n -= d; return *this; }

bool rational::operator==(const rational& q) const{
	if(is_nan() || q.is_nan()){ return false; }
	return ((q.n == n) && (q.d == d));
}
bool rational::operator!=(const rational& q) const{
	if(is_nan() || q.is_nan()){ return false; }
	return ((q.n != n) || (q.d != d));
}
bool rational::operator< (const rational& q) const{
	if(is_nan() || q.is_nan()){ return false; }

	// Determine relative order by expanding each value to its simple continued
	// fraction representation using the Euclidian GCD algorithm.
	struct { int n, d, q, r; }
		ts = { this->n, this->d, this->n / this->d, this->n % this->d },
		rs = { q.n, q.d, q.n / q.d, q.n % q.d };
	unsigned int reverse = 0;

	// Normalize negative moduli by repeatedly adding the (positive) denominator
	// and decrementing the quotient.  Later cycles should have all positive
	// values, so this only has to be done for the first cycle.  (The rules of
	// C++ require a nonnegative quotient & remainder for a nonnegative dividend
	// & positive divisor.)
	while(ts.r < 0){ ts.r += ts.d; --ts.q; }
	while(rs.r < 0){ rs.r += rs.d; --rs.q; }

	// Loop through and compare each variable's continued-fraction components
	while(true){
		// The quotients of the current cycle are the continued-fraction
		// components.  Comparing two c.f. is comparing their sequences,
		// stopping at the first difference.
		if ( ts.q != rs.q ){
			// Since reciprocation changes the relative order of two variables,
			// and c.f. use reciprocals, the less/greater-than test reverses
			// after each index.  (Start w/ non-reversed @ whole-number place.)
			return reverse ? ts.q > rs.q : ts.q < rs.q;
		}

		// Prepare the next cycle
		reverse ^= 1;

		if((0 == ts.r) || (0 == rs.r)){
			// At least one variable's c.f. expansion has ended
			break;
		}

		ts.n = ts.d;         ts.d = ts.r;
		ts.q = ts.n / ts.d;  ts.r = ts.n % ts.d;
		rs.n = rs.d;         rs.d = rs.r;
		rs.q = rs.n / rs.d;  rs.r = rs.n % rs.d;
	}

	// Compare infinity-valued components for otherwise equal sequences
	if(ts.r == rs.r){
		// Both remainders are zero, so the next (and subsequent) c.f.
		// components for both sequences are infinity.  Therefore, the sequences
		// and their corresponding values are equal.
		return false;
	}else{
		// Exactly one of the remainders is zero, so all following c.f.
		// components of that variable are infinity, while the other variable
		// has a finite next c.f. component.  So that other variable has the
		// lesser value (modulo the reversal flag!).
		return (ts.r != 0) != (reverse != 0);
	}
}
bool rational::operator<=(const rational& q) const{
	if(is_nan() || q.is_nan()){ return false; }
	return ((*this == q) || (*this < q));
}
bool rational::operator> (const rational& q) const{
	if(is_nan() || q.is_nan()){ return false; }
	return ((*this != q) && !(*this < q));
}
bool rational::operator>=(const rational& q) const{
	if(is_nan() || q.is_nan()){ return false; }
	return ((*this == q) || !(*this < q));
}

bool rational::operator==(int numerator) const{
	if(is_nan()){ return false; }
	return ((numerator == n) && (1 == d));
}
bool rational::operator!=(int numerator) const{
	if(is_nan()){ return false; }
	return ((numerator != n) || (1 != d));
}
bool rational::operator< (int numerator) const{
	if(is_nan()){ return false; }

	// Break value into mixed-fraction form, w/ always-nonnegative remainder
	int q = this->n / this->d, r = this->n % this->d;
	while(r < 0){ r += this->d; --q; }

	// Compare with just the quotient, since the remainder always bumps the
	// value up.  [Since q = floor(n/d), and if n/d < i then q < i, if n/d == i
	// then q == i, if n/d == i + r/d then q == i, and if n/d >= i + 1 then
	// q >= i + 1 > i; therefore n/d < i iff q < i.]
	return q < numerator;
}
bool rational::operator<=(int numerator) const{
	if(is_nan()){ return false; }
	return ((*this == numerator) || (*this < numerator));
}
bool rational::operator> (int numerator) const{
	if(is_nan()){ return false; }
	return ((*this != numerator) && !(*this < numerator));
}
bool rational::operator>=(int numerator) const{
	if(is_nan()){ return false; }
	return ((*this == numerator) || !(*this < numerator));
}

float rational::float_value() const{ return (float)n/(float)d; }
double rational::double_value() const{ return (double)n/(double)d; }

unsigned int rational::gcd(unsigned int u, unsigned int v){ // from wikipedia
	int shift;

	// GCD(0,x) := x
	if (u == 0 || v == 0)
	return u | v;

	// Let shift := lg K, where K is the greatest power of 2
	// dividing both u and v.
	for(shift = 0; ((u | v) & 1) == 0; ++shift){
		u >>= 1;
		v >>= 1;
	}

	while((u & 1) == 0){ u >>= 1; }

	// From here on, u is always odd.
	do{
		while((v & 1) == 0){ v >>= 1; }

		// Now u and v are both odd, so diff(u, v) is even.
		// Let u = min(u, v), v = diff(u, v)/2.
		if(u < v){
			v -= u;
		}else{
			unsigned int diff = u - v;
			u = v;
			v = diff;
		}
		v >>= 1;
	}while(v != 0);

	return u << shift;
}
unsigned int rational::lcm(unsigned int u, unsigned int v){
	unsigned int g = rational::gcd(u, v);
	if(0 == g){ return 0; }
	return (u/g)*v;
}
int rational::isqrt(int u){
	if(u <= 0){ return 0; }
	int a = 2;
	int b = u/a;
	while(((a-b) > 1) || ((b-a) > 1)){
		a = (a+b)/2;
		b = u/a;
	}
	return (a<b)?a:b;
}
double rational::double_to_rational(double x, int &numerator, int &denominator){
	/*
	** find rational approximation to given real number
	** David Eppstein / UC Irvine / 8 Aug 1993
	**
	** With corrections from Arno Formella, May 2008
	**
	** usage: a.out r d
	**	 r is real number to approx
	**	 d is the maximum denominator allowed
	**
	** based on the theory of continued fractions
	** if x = a1 + 1/(a2 + 1/(a3 + 1/(a4 + ...)))
	** then best approximation is found by truncating this series
	** (with some adjustments in the last term).
	**
	** Note the fraction can be recovered as the first column of the matrix
	**	( a1 1 ) ( a2 1 ) ( a3 1 ) ...
	**	( 1	 0 ) ( 1  0 ) ( 1  0 )
	** Instead of keeping the sequence of continued fraction terms,
	** we just keep the last partial product of these matrices.
	*/
	int m[2][2];
	int sign = 0;
	if(x < 0){ x = -x; sign = 1; }
	double startx = x;
	const int maxden = 10000;
	int ai;

	// initialize matrix
	m[0][0] = m[1][1] = 1;
	m[0][1] = m[1][0] = 0;

	// loop finding terms until denom gets too big
	while(m[1][0] *  ( ai = (int)x ) + m[1][1] <= maxden){
		int t;
		t = m[0][0] * ai + m[0][1];
		m[0][1] = m[0][0];
		m[0][0] = t;
		t = m[1][0] * ai + m[1][1];
		m[1][1] = m[1][0];
		m[1][0] = t;
		if(x == (double)ai) break;	 // AF: division by zero
		x = 1/(x - (double) ai);
		if(x > (double)0x7FFFFFFF) break;	 // AF: representation failure
	}

	// now remaining x is between 0 and 1/ai
	// approx as either 0 or 1/m where m is max that will fit in maxden
	// first try zero
	numerator = m[0][0];
	denominator = m[1][0];
	double error = startx - (double)numerator/(double)denominator;
	
	// now try other possibility
	ai = (maxden - m[1][1]) / m[1][0];
	m[0][0] = m[0][0] * ai + m[0][1];
	m[1][0] = m[1][0] * ai + m[1][1];
	int n2 = m[0][0], d2 = m[1][0];
	double e2 = startx - ((double) m[0][0] / (double) m[1][0]);
	if(fabs(e2) <fabs(error)){
		numerator = n2;
		denominator = d2;
		error = e2;
	}
	if(sign){ numerator = -numerator; }
	return error;
}

rational operator+(const rational& a, const rational &b){
	rational ret(a);
	ret += b;
	return ret;
}
rational operator-(const rational& a, const rational &b){
	rational ret(a);
	ret -= b;
	return ret;
}
rational operator*(const rational& a, const rational &b){
	rational ret(a);
	ret *= b;
	return ret;
}
rational operator/(const rational& a, const rational &b){
	rational ret(a);
	ret /= b;
	return ret;
}
rational operator-(const rational& q){
	return rational(-q.num(), q.den());
}

rational abs(const rational& q){
	return rational(abs(q.num()), q.den());
}
int floor(const rational& q){
	if(q.is_nan()){ return 0; }
	int n = q.num(), d = q.den();
	if(n > 0){
		return n/d;
	}else if(0 == n){
		return 0;
	}else{
		return -((-n)/d)-1;
	}
}
int ceil(const rational& q){
	if(q.is_nan()){ return 0; }
	int n = q.num(), d = q.den();
	if(0 == n){
		return 0;
	}else if(0 == n%d){
		return n/d;
	}else if(n > 0){
		return (n/d)+1;
	}else{
		return -((-n)/d);
	}
}
int round(const rational& q){
	int fq = floor(q);
	rational r = q - fq;
	if(r > rational(1,2)){
		return fq+1;
	}else{
		return fq;
	}
}


std::istream& operator>>(std::istream& is, rational& q){
	int n = 0, d = 1;
	char c = 0;
	std::istream::fmtflags oldflags = is.flags();

	is >> n;
	c = is.get();

	if(c != '/'){ is.clear(std::istream::badbit); }

#if !defined(__GNUC__) || (defined(__GNUC__) && (__GNUC__ >= 3)) || defined __SGI_STL_PORT
	is >> std::noskipws;
#else
	is.unsetf(ios::skipws); // compiles, but seems to have no effect.
#endif
	is >> d;

	if(is){ q = rational(n, d); }

	is.flags(oldflags);
	return is;
}
std::ostream& operator<<(std::ostream& os, const rational& q){
	if(q.is_nan()){
		os << "(nan)";
	}else if(1 == q.den()){
		os << q.num();
	}else{
		os << q.num() << '/' << q.den();
	}
	return os;
}


#ifdef USING_NUMERIC_TYPE_TRAITS
template <>
float ScalarTraits<rational>::numeric_value<float>(const value_type &v){ return v.float_value(); }
#endif

