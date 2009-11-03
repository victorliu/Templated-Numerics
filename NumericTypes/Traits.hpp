#ifndef _NUMERIC_TYPE_TRAITS_H_
#define _NUMERIC_TYPE_TRAITS_H_

#include <limits>
#include <complex>

#define USING_NUMERIC_TYPE_TRAITS

template <typename NumericType>
class ScalarTraits{
public:
	typedef NumericType value_type;
	
	//template <typename FloatType>
	//static FloatType numeric_value(const value_type &v);
	
	// Must have these:
	//static const value_type EPSILON;
	//static const value_type MAX_VALUE;
	//static const value_type MIN_VALUE;
};

template <typename NumericType>
class ComplexTraits{
public:
	typedef NumericType value_type;
	// Must define this:
	//typedef NumericType real_t; // or whatever
	
	virtual value_type conj(const value_type &a) = 0;
	virtual value_type real(const value_type &a) = 0;
	virtual value_type imag(const value_type &a) = 0;
	
	// Must have these:
	//static const value_type ZERO;
	//static const value_type ONE;
};

template <typename NumericType>
class FieldTraits{
public:
	typedef NumericType value_type;
	
	// Must define this:
	//static value_type Solve(const value_type &A, const value_type &b, value_type &x); // solves Ax=b
	
	// Must have these:
	//static const value_type ZERO;
	//static const value_type ONE;
};

//////////////////////////////// Specializations

template <>
class ScalarTraits<double>{
public:
	typedef double value_type;
	
	template <typename FloatType>
	static FloatType numeric_value(const value_type &v){ return (FloatType)v; }
	
	static const value_type epsilon(){ return std::numeric_limits<double>::epsilon(); }
	static const value_type max_value(){ return std::numeric_limits<double>::max(); }
	static const value_type min_value(){ return std::numeric_limits<double>::min(); }
};

template <>
class ScalarTraits<float>{
public:
	typedef float value_type;
	
	template <typename FloatType>
	static FloatType numeric_value(const value_type &v){ return (FloatType)v; }
	
	static const value_type epsilon(){ return std::numeric_limits<float>::epsilon(); }
	static const value_type max_value(){ return std::numeric_limits<float>::max(); }
	static const value_type min_value(){ return std::numeric_limits<float>::min(); }
};


template <>
class ComplexTraits<std::complex<double> >{
public:
	typedef std::complex<double> value_type;
	typedef double real_t;
	
	value_type conj(const value_type &a){ return std::conj(a); }
	value_type real(const value_type &a){ return a.real(); }
	value_type imag(const value_type &a){ return a.imag(); }
	
	static const value_type zero(){ return value_type(0); }
	static const value_type one(){ return value_type(1); }
};


template <>
class ComplexTraits<std::complex<float> >{
public:
	typedef std::complex<float> value_type;
	typedef double real_t;
	
	value_type conj(const value_type &a){ return std::conj(a); }
	value_type real(const value_type &a){ return a.real(); }
	value_type imag(const value_type &a){ return a.imag(); }
	
	static const value_type zero(){ return value_type(0); }
	static const value_type one(){ return value_type(1); }
};



template <>
class FieldTraits<double>{
public:
	typedef double value_type;
	
	static value_type Solve(const value_type &A, const value_type &b, value_type &x){
		x = b/A;
	}
	
	static const value_type zero(){ return value_type(0); }
	static const value_type one(){ return value_type(1); }
};

template <>
class FieldTraits<float>{
public:
	typedef float value_type;
	
	static value_type Solve(const value_type &A, const value_type &b, value_type &x){
		x = b/A;
	}
	
	static const value_type zero(){ return value_type(0); }
	static const value_type one(){ return value_type(1); }
};

template <>
class FieldTraits<std::complex<double> >{
public:
	typedef std::complex<double> value_type;
	
	static value_type Solve(const value_type &A, const value_type &b, value_type &x){
		x = b/A;
	}
	
	static const value_type zero(){ return value_type(0); }
	static const value_type one(){ return value_type(1); }
};

template <>
class FieldTraits<std::complex<float> >{
public:
	typedef std::complex<float> value_type;
	
	static value_type Solve(const value_type &A, const value_type &b, value_type &x){
		x = b/A;
	}
	
	static const value_type zero(){ return value_type(0); }
	static const value_type one(){ return value_type(1); }
};

#endif // _NUMERIC_TYPE_TRAITS_H_

