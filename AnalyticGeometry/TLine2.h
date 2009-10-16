#ifndef _TLINE2_H_
#define _TLINE2_H_

#include "TPt2.h"
#include "TVec2.h"

template <class NumericType>
class TLine2{
public:
	typedef NumericType value_type;
	typedef TPt2<value_type> Pt2;
	typedef TVec2<value_type> Vec2;

	Pt2 p;
	Vec2 v;
	
	TLine2(const Pt2 &A, const Pt2 &B):p(A),v(B-A){}
	TLine2(const Pt2 &A, const Vec2 &AB):p(A),v(AB){}
	TLine2(const TLine2 &l):p(l.p),v(l.v){}
	
	// returns 1 if a is on left of the line (stand at p, look towards v)
	// returns -1 if a is on right, 0 if exactly on line
	int SideOf(const Pt2 &a) const{
		NumericType c(Vec2::Cross(v, a-p));
		if(0 == c){ return 0; }
		else if(c > 0){ return 1; }
		else{ return -1; }
	}
	Pt2 operator[](const NumericType &t) const{
		return p+t*v;
	}
	// Inverse of operator[]; projects a on line and returns the parametric coordinate
	NumericType Projection(const Pt2 &a){
		return Vec2::Dot(v,a-p)/v.LengthSq();
	}
};

#endif // _TLINE2_H_
