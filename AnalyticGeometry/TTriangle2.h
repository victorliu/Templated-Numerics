#ifndef _TTRIANGLE2_H_
#define _TTRIANGLE2_H_

#include "TPt2.h"
#include "TVec2.h"

template <class NumericType>
class TTriangle2{
public:
	typedef NumericType value_type;
	typedef TPt2<value_type> Pt2;
	typedef TVec2<value_type> Vec2;

	Pt2 base;
	Vec2 u, v;
	
	// Orientation is by no means guaranteed just because the constructors enforce positive orientation
	// If you mess with the members by yourself, that's your own problem; you call Orient() yourself.
	TTriangle2(const Pt2 &Base, const Vec2 &edge1, const Vec2 &edge2):base(Base),u(edge1),v(edge2){
		Orient();
	}
	TTriangle2(const Pt2 &a, const Pt2 &b, const Pt2 &c):base(a),u(b-a),v(c-a){
		Orient();
	}
	TTriangle2(const TTriangle2 &p):base(p.base),u(p.u),v(p.v){}
	
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
		which = which % 3;
		switch(which){
		case 0:
			return base;
		case 1:
			return base+u;
		default:
			return base+(u+v);
		}
	}
	value_type Area() const{
		return (value_type(1)/value_type(2)) * abs(Vec2::Cross(u,v));
	}
	// The vectors will sum to 1, but orientation is arbitrary
	Vec2 Side(size_t which) const{
		which = which % 3;
		switch(which){
		case 0:
			return u;
		case 1:
			return (v-u);
		default:
			return (-v);
		}
	}
	void Orient(){ // ensures positive orientation
		if(Vec2::Cross(u,v) < 0){
			std::swap(u,v);
		}
	}
	value_type LongestSideLength() const{
		value_type len[3];
		len[0] = u.Length();
		len[1] = v.Length();
		len[2] = (v-u).Length();
		if(len[0] > len[1]){
			if(len[0] > len[2]){
				return len[0];
			}else{
				return len[2];
			}
		}else{
			if(len[1] > len[2]){
				return len[1];
			}else{
				return len[2];
			}
		}
	}
	value_type ShortestSideLength() const{
		value_type len[3];
		len[0] = u.Length();
		len[1] = v.Length();
		len[2] = (v-u).Length();
		if(len[0] < len[1]){
			if(len[0] < len[2]){
				return len[0];
			}else{
				return len[2];
			}
		}else{
			if(len[1] < len[2]){
				return len[1];
			}else{
				return len[2];
			}
		}
	}
};

#endif // _TTRIANGLE2_H_
