#ifndef _TPARALLELOGRAM2_H_
#define _TPARALLELOGRAM2_H_

#include "TPt2.h"
#include "TVec2.h"

template <class NumericType>
class TParallelogram2{
public:
	typedef NumericType value_type;
	typedef TPt2<value_type> Pt2;
	typedef TVec2<value_type> Vec2;

	Pt2 base;
	Vec2 u, v;
	
	TParallelogram2(const Pt2 &Base, const Vec2 &edge1, const Vec2 &edge2):base(Base),u(edge1),v(edge2){}
	TParallelogram2(const TParallelogram2 &p):base(p.base),u(p.u),v(p.v){}
	
	bool Contains(const Pt2 &p) const{
		// The parallel sides define hyperplanes
		// A point is inside if it is not on the same side of both hyperplanes
		{
			Vec2 n(u.Rotate90CCW());
			value_type a = Vec2::Dot(p-base, n);
			value_type b = Vec2::Dot(p-(base+v), n);
			if((a > 0 && b > 0) || (a < 0 && b < 0)){ return false; }
		}
		{
			Vec2 n(v.Rotate90CCW());
			value_type a = Vec2::Dot(p-base, n);
			value_type b = Vec2::Dot(p-(base+u), n);
			if((a > 0 && b > 0) || (a < 0 && b < 0)){ return false; }
		}
		return true;
	}
	Pt2 operator[](size_t which) const{
		switch(which){
		case 1:
			return base+u+v;
		case 2:
			return base+u;
		case 3:
			return base+v;
		default:
			return base;
		}
	}
	value_type Area() const{
		return abs(Vec2::Cross(u,v));
	}
};

#endif // _TPARALLELOGRAM2_H_
