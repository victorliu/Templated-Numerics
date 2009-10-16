#ifndef _TAABB2_H_
#define _TAABB2_H_

#include "TPt2.h"
#include "TVec2.h"

template <class NumericType>
class TAABB2{
public:
	typedef NumericType value_type;
	typedef TPt2<value_type> Pt2;
	typedef TVec2<value_type> Vec2;

	Pt2 p0, p1;
	
	TAABB2(const Pt2 &P0, const Pt2 &P1):p0(P0),p1(P1){
		// check to make sure p0.r[i] < p1.r[i]
		if(p0.r[0] > p1.r[0]){ std::swap(p0.r[0], p1.r[0]); }
		if(p0.r[1] > p1.r[1]){ std::swap(p0.r[1], p1.r[1]); }
	}
	TAABB2(const TAABB2 &b):p0(b.p0),p1(b.p1){}
	
	bool Contains(const Pt2 &p) const{
		return p0.r[0] <= p.r[0] && p.r[0] <= p1.r[0] && p0.r[1] <= p.r[1] && p.r[1] <= p1.r[1];
	}
};

#endif // _TAABB2_H_
