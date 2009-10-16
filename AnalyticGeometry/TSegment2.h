#ifndef _TSEGMENT2_H_
#define _TSEGMENT2_H_

#include "TPt2.h"
#include "TVec2.h"

template <class NumericType>
class TSegment2{
public:
	typedef NumericType value_type;
	typedef TPt2<value_type> Pt2;
	typedef TVec2<value_type> Vec2;

	Pt2 a,b;
	
	TSegment2(const Pt2 &A, const Pt2 &B):a(A),b(B){}
	TSegment2(const Pt2 &A, const Vec2 &AB):a(A),b(A+AB){}
	TSegment2(const TSegment2 &s):a(s.a),b(s.b){}
	
	const Pt2& operator[](size_t which) const{
		return ((which & 1) == 0) ? a : b;
	}
	Pt2& operator[](size_t which){
		return ((which & 1) == 0) ? a : b;
	}
	value_type Length() const{
		return (b-a).Length();
	}
};

#endif // _TSEGMENT2_H_
