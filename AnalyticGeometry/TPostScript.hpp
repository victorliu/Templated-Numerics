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
	static const NumericType cx(612/2), cy(792/2);
	//static const NumericType cx(0), cy(0);
	return TPt2<NumericType>(
		cx+TransformLinear<NumericType,InputType>(p.r[0]),
		cy+TransformLinear<NumericType,InputType>(p.r[1])
		);
}
template <typename NumericType, typename InputType>
TVec2<NumericType> TransformLinear(const TVec2<InputType> &v){
	return TVec2<NumericType>(
		TransformLinear<NumericType,InputType>(v.v[0]),
		TransformLinear<NumericType,InputType>(v.v[1])
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


template <typename InputType>
void Initialize(std::ostream &os = std::cout){
	os << "/Times-Roman findfont" << std::endl;
	os << "6 scalefont" << std::endl;
	os << "setfont" << std::endl;
}
template <typename InputType>
void InitializeEPS(std::ostream &os = std::cout){
	os << "%!PS-Adobe EPSF-3.0" << std::endl;
	os << "%%Creator: TPostScript.hpp" << std::endl;
	Initialize<InputType>(os);
}
template <typename InputType>
void Close(std::ostream &os = std::cout){
	os << "showpage" << std::endl;
}
template <typename InputType>
void CloseEPS(const TPt2<InputType> &minpt, const TPt2<InputType> &maxpt, std::ostream &os = std::cout){
	os << "%%BoundingBox: " << minpt.r[0] << " " << minpt.r[1] << " " << maxpt.r[0] << " " << maxpt.r[1] << std::endl;
	os << "%%Origin: " << minpt.r[0] << " " <<minpt.r[1] << std::endl;
	os << "%%LanguageLevel: 2" << std::endl;
	os << "%%Pages: 1" << std::endl;
	os << "%%Page: 1 1" << std::endl;
	os << "%%EOF" << std::endl;
}

template <typename InputType>
void SetColor(const InputType &red, const InputType &green, const InputType &blue, std::ostream &os = std::cout){
	os << red << " " << green << " " << blue << " setrgbcolor" << std::endl;
}

#include "TPt2.h"
#include "TVec2.h"

template <typename InputType>
void DrawText(const TPt2<InputType> &lower_left_point, std::string str, std::ostream &os = std::cout){
	TPt2<float> p1(TransformAffine<float,InputType>(lower_left_point));
	os << "newpath" << std::endl;
	os << p1.r[0] << " " << p1.r[1] << " moveto" << std::endl;
	// All parentheses must be escaped:
	size_t pos = 0;
	while((pos = str.find("(", pos)) != std::string::npos){
		str.replace(pos, 1, "\\(");
		pos += 2;
	}
	pos = 0;
	while((pos = str.find(")", pos)) != std::string::npos){
		str.replace(pos, 1, "\\)");
		pos += 2;
	}

	os << "(" << str << ") show" << std::endl;
}

template <typename InputType>
void DrawPoint(const TPt2<InputType> &p, std::ostream &os = std::cout){
	static const float rad = 1.5;
	/*
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
	*/
	TPt2<float> c(TransformAffine<float,InputType>(p));
	os << "newpath" << std::endl;
	os << c.r[0] << ' ' << c.r[1] << ' ' << rad << " 0 360 arc closepath" << std::endl;
	os << "fill" << std::endl;
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
void DrawArrow(const TPt2<InputType> &base, const TVec2<InputType> &vec, std::ostream &os = std::cout){
	TPt2<float> b(TransformAffine<float,InputType>(base));
	TVec2<float> v(TransformLinear<float,InputType>(vec));
	TPt2<float> t(TransformAffine<float,InputType>(base+vec));
	TVec2<float> n(0.05f*TVec2<float>::Rotate90CCW(v)); // normal to vector
	TPt2<float> s0(b + 0.9f*v + n);
	TPt2<float> s1(b + 0.9f*v - n);
	
	os << "newpath" << std::endl;
	os << b.r[0] << ' ' << b.r[1] << " moveto" << std::endl;
	os << t.r[0] << ' ' << t.r[1] << " lineto" << std::endl;
	os << "stroke" << std::endl;
	os << "newpath" << std::endl;
	os << s0.r[0] << ' ' << s0.r[1] << " moveto" << std::endl;
	os << t.r[0] << ' ' << t.r[1] << " lineto" << std::endl;
	os << s1.r[0] << ' ' << s1.r[1] << " lineto" << std::endl;
	os << "stroke" << std::endl;
}
template <typename InputType>
void DrawCircle(const TPt2<InputType> &center, const InputType &radius, std::ostream &os = std::cout){
	TPt2<float> c = TransformAffine<float,InputType>(center);
	os << c.r[0] << ' ' << c.r[1] << ' ' << TransformLinear<float,InputType>(radius) << " 0 360 arc closepath" << std::endl;
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

