#ifndef _TCIRCLE2_H_
#define _TCIRCLE2_H_

#include "TPt2.h"

template <class NumericType>
class TCircle2{
public:
	typedef NumericType value_type;
	typedef TPt2<value_type> Pt2;

	Pt2 center;
	value_type radius;

	TCircle2(const Pt2 &c, const value_type &r):center(c),radius(r){}
	TCircle2(const TCircle2 &c):center(c.center),radius(c.radius){}

	bool Contains(const Pt2 &p) const{
		return (p-center).LengthSq() < (radius*radius);
	}
	value_type Area() const{
		return value_type(M_PI) * radius*radius;
	}
};

#endif // _TCIRCLE2_H_
