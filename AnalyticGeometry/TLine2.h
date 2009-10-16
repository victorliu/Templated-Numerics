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
};

#endif // _TLINE2_H_
