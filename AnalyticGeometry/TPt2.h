#ifndef _TPT2_H_
#define _TPT2_H_

#include "TVec2.h"

template <class RealType>
struct TPt2{
	typedef RealType value_type;
	typedef RealType real_type;
	typedef TPt2<RealType> Pt2;
	real_type r[2];
	
	static const TPt2 Origin;
	
	TPt2(){ r[0] = 0; r[1] = 0; }
	TPt2(const real_type &x, const real_type &y){
		r[0] = x;
		r[1] = y;
	}
	TPt2(const TPt2 &p){ r[0] = p.r[0]; r[1] = p.r[1]; }
	
	const real_type& operator()(size_t idx) const{ return r[idx]; }
	real_type& operator()(size_t idx){ return r[idx]; }
	const real_type& x() const{ return r[0]; }
	const real_type& y() const{ return r[1]; }
	real_type& x(){ return r[0]; }
	real_type& y(){ return r[1]; }
	
	TPt2<real_type>& operator += (const TVec2<real_type> &a){ r[0] += a.v[0]; r[1] += a.v[1]; }
	TPt2<real_type>& operator -= (const TVec2<real_type> &a){ r[0] -= a.v[0]; r[1] -= a.v[1]; }
};
template <class RealType>
const TPt2<RealType> TPt2<RealType>::Origin(0,0);

template <class RealType>
TVec2<RealType> operator - (const TPt2<RealType> &a, const TPt2<RealType> &b){
	return TVec2<RealType>(a.r[0]-b.r[0], a.r[1]-b.r[1]);
}

template <class RealType>
TPt2<RealType> operator + (const TPt2<RealType> &a, const TVec2<RealType> &b){
	return TPt2<RealType>(a.r[0]+b.v[0], a.r[1]+b.v[1]);
}

template <class RealType>
TPt2<RealType> operator - (const TPt2<RealType> &a, const TVec2<RealType> &b){
	return TPt2<RealType>(a.r[0]-b.v[0], a.r[1]-b.v[1]);
}
template <class RealType>
std::ostream& operator<<(std::ostream &os, TPt2<RealType> &p){
	os << '{' << p.r[0] << ", " << p.r[1] << ')';
}

template <class RealType>
struct TPt2LexicalLess{ bool operator()(const TPt2<RealType> &p1, const TPt2<RealType> &p2) const{
	if(p1.r[0] < p2.r[0]){ return true; }
	else if(p1.r[0] == p2.r[0]){
		return p1.r[1] < p2.r[1];
	}else{
		return false;
	}
}};

#endif // _TPT2_H_
