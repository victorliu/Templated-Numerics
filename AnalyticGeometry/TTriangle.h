#ifndef _TTRIANGLE_H_
#define _TTRIANGLE_H_

#include "TPt2.h"
#include "TVec2.h"

template <class NumericType>
class TTriangle{
public:
	typedef NumericType value_type;
	typedef TPt2<value_type> Pt2;
	typedef TVec2<value_type> Vec2;

	Pt2 base;
	Vec2 u, v;
	
	TTriangle(const Pt2 &Base, const Vec2 &edge1, const Vec2 &edge2):base(Base),u(edge1),v(edge2){}
	TTriangle(const Pt2 &a, const Pt2 &b, const Pt2 &c):base(a),u(b-a),v(c-a){}
	TTriangle(const TTriangle &p):base(p.base),u(p.u),v(p.v){}
	
	bool Contains(const Pt2 &p) const{
		// We transform into the unit orthogonal basis and check the coords
		// [s;t] = [u v]^{-1} (p-base)
		// s >= 0, t >= 0, s+t <= 1
		// We compute s and t before dividing through by the det
		Vec2 d(p-base);
		value_type s = v[1]*d[0] - v[0]*d[1];
		value_type t = u[0]*d[0] - u[1]*d[1];
		if(s < 0 || t < 0){ return false; }
		value_type det = u[0]*v[1] - u[1]*v[0];
		return (s + t <= det);
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
		return (value_type(1)/value_type(2)) * abs(Vec2::Cross(u,v));
	}
};

#endif // _TTRIANGLE_H_
