#ifndef _INTEGRATION_HPP_
#define _INTEGRATION_HPP_

namespace Integration{

enum Method{
	MIDPOINT,
	TRAPEZOIDAL,
	SIMPSON,
	GAUSS_KRONROD
};

template <typename NumericType>
struct Parameters{
	size_t max_evaluations;
	NumericType max_relative_error;
	NumericType max_absolute_error;
};

template <typename NumericType>
class Function1{
public:
	virtual NumericType operator(const NumericType &x) = 0;
};

// Forward declarations for all 1D integration methods
template <typename NumericType>
NumericType Integrate_Midpoint(const Function1 &f, const Parameters &params);
template <typename NumericType>
NumericType Integrate_Trapezoidal(const Function1 &f, const Parameters &params);
template <typename NumericType>
NumericType Integrate_Simpson(const Function1 &f, const Parameters &params);
template <typename NumericType>
NumericType Integrate_GaussKronrod(const Function1 &f, const Parameters &params);

template <typename NumericType>
NumericType Integrate(const Function1 &f, const Parameters &params, Method method = MIDPOINT){
	switch(method){
	case MIDPOINT:
		return Integrate_Midpoint(f, params);
	case TRAPEZOIDAL:
		return Integrate_Trapezoidal(f, params);
	case SIMPSON:
		return Integrate_Simpson(f, params);
	case GAUSS_KRONROD:
		return Integrate_GaussKronrod(f, params);
	default: // MIDPOINT
		return Integrate_Midpoint(f, params);
	}
}

};

#endif // _INTEGRATION_HPP_
