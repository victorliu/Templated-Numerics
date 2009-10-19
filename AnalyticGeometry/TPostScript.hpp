#ifndef _TPOSTSCRIPT_HPP_
#define _TPOSTSCRIPT_HPP_

// Preprocessor switches:
//   TPOSTSCRIPT_LIGHT - only define and use a minimum number of primitives

namespace TPostScript{

// Transformation to page coordinates

template <typename NumericType>
NumericType TransformLinearNumeric(NumericType x){
	static const NumericType scale(40);
	return scale*x;
}
template <typename NumericType, typename InputType>
NumericType TransformLinear(InputType x){
	return TransformLinearNumeric(x.numeric_value());
}
template <typename NumericType, typename InputType>
TPt2<NumericType> TransformAffine(const TPt2<InputType> &p){
	//static const NumericType cx(612/2), cy(792/2);
	static const NumericType cx(0), cy(0);
	return TPt2<NumericType>(
		cx+TransformLinear<NumericType,InputType>(p.r[0]),
		cy+TransformLinear<NumericType,InputType>(p.r[1])
		);
}
template <typename NumericType, typename InputType>
TVec2<NumericType> TransformLinear(const TVec2<InputType> &p){
	return TVec2<NumericType>(
		TransformLinear<NumericType,InputType>(p.r[0]),
		TransformLinear<NumericType,InputType>(p.r[1])
		);
}

// Specializations for primitive input types
template <>
double TransformLinear<double,double>(double x){
	return TransformLinearNumeric(double(x));
}
template <>
float TransformLinear<float,float>(float x){
	return TransformLinearNumeric(float(x));
}
template <>
double TransformLinear<double,float>(float x){
	return TransformLinearNumeric(double(x));
}
template <>
float TransformLinear<float,double>(double x){
	return TransformLinearNumeric(float(x));
}


#include "TPt2.h"
#include "TVec2.h"

template <typename InputType>
void DrawPoint(const TPt2<InputType> &p, std::ostream &os = std::cout){
	static const float rad = 1;
	os << "newpath" << std::endl;
	TPt2<float> p0, p1;
	p0 = TransformAffine<float,InputType>(p); p1 = p0;
	p0.r[0] -= rad; p1.r[0] += rad;
	os << p0.r[0] << ' ' << p0.r[1] << " moveto" << std::endl;
	os << p1.r[0] << ' ' << p1.r[1] << " lineto" << std::endl;
	os << "stroke" << std::endl;
	
	os << "newpath" << std::endl;
	p0 = TransformAffine<float,InputType>(p); p1 = p0;
	p0.r[1] -= rad; p1.r[1] += rad;
	os << p0.r[0] << ' ' << p0.r[1] << " moveto" << std::endl;
	os << p1.r[0] << ' ' << p1.r[1] << " lineto" << std::endl;
	os << "stroke" << std::endl;
}
template <typename InputType>
void DrawSegment(const TPt2<InputType> &a, const TPt2<InputType> &b, std::ostream &os = std::cout){
	TPt2<float> p0(TransformAffine<float,InputType>(a));
	TPt2<float> p1(TransformAffine<float,InputType>(b));
	
	os << "newpath" << std::endl;
	os << p0.r[0] << ' ' << p0.r[1] << " moveto" << std::endl;
	os << p1.r[0] << ' ' << p1.r[1] << " lineto" << std::endl;
	os << "stroke" << std::endl;
}
template <typename InputType>
void DrawCircle(const TPt2<InputType> &center, const InputType &radius, std::ostream &os = std::cout){
	TPt2<float> c = TransformAffine<float,InputType>(center);
	os << c.r[0] << ' ' << c.r[1] << ' ' << TransformLinear(radius) << " 0 360 arc closepath" << std::endl;
	os << "stroke" << std::endl;
}

#ifndef TPOSTSCRIPT_LIGHT

#include "TSegment2.h"
#include "TCircle2.h"

template <typename InputType>
void DrawCircle(const TCircle<InputType> &c, std::ostream &os = std::cout){
	DrawCircle(c.center, c.radius, os);
}

template <typename InputType>
void DrawSegment(const TSegment<InputType> &s, std::ostream &os = std::cout){
	DrawSegment(s[0], s[1], os);
}

#endif // TPOSTSCRIPT_LIGHT

}; // namespace TPostScript

#endif // _TPOSTSCRIPT_HPP_

