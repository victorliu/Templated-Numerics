#ifndef _TVEC2_H_
#define _TVEC2_H_

#include <iostream>

template <class RealType>
struct TVec2{
	typedef RealType value_type;
	typedef RealType real_type;
	real_type v[2];
	
	static const TVec2 Zero;
	
	TVec2(){ v[0] = 0; v[1] = 0; }
	TVec2(const real_type &x, const real_type &y){
		v[0] = x;
		v[1] = y;
	}
	TVec2(const TVec2 &p){ v[0] = p.v[0]; v[1] = p.v[1]; }
	TVec2& operator=(const TVec2 &p){ v[0] = p.v[0]; v[1] = p.v[1]; return *this; }
	
	const real_type& operator[](size_t idx) const{ return v[idx]; }
	real_type& operator[](size_t idx){ return v[idx]; }
	const real_type& x() const{ return v[0]; }
	const real_type& y() const{ return v[1]; }
	real_type& x(){ return v[0]; }
	real_type& y(){ return v[1]; }	/// Normalizes the vector. If zero vector, produces an infinite vector.

	void Normalize(){
		real_type L = Length();
		v[0] /= L;
		v[1] /= L;
	}

	/// Returns the L^2 norm of the vector.
	real_type Length() const{
		if(abs(v[1]) > abs(v[0])){
			return v[1]*sqrt(real_type(1) + v[0]/v[1]);
		}else{
			return v[0]*sqrt(real_type(1) + v[1]/v[0]);
		}
	}
	/// Returns the square of the L^2 norm of the vector.
	real_type LengthSq() const{
		return v[0]*v[0] + v[1]*v[1];
	}

	/// Performs the dot (inner) product on two vectors (Euclidean R^2 metric).
	static real_type Dot(const TVec2<real_type> &a, const TVec2<real_type> &b){
		return a.v[0]*b.v[0] + a.v[1]*b.v[1];
	}
	
	/// Performs the scalar product (length of cross product).
	static real_type Cross(const TVec2<real_type> &a, const TVec2<real_type> &b){
		return a.v[0]*b.v[1] - a.v[1]*b.v[0];
	}
	static TVec2<real_type> Rotate90CCW(const TVec2<real_type> &v){
		return TVec2<real_type>(-v[1], v[0]);
	}

	TVec2<real_type>& operator *= (const real_type &C){ v[0] *= C; v[1] *= C; return *this; }
	TVec2<real_type>& operator /= (const real_type &C){ v[0] /= C; v[1] /= C; return *this; }
	TVec2<real_type>& operator += (const TVec2<real_type> &a){ v[0] += a.v[0]; v[1] += a.v[1]; return *this; }
	TVec2<real_type>& operator -= (const TVec2<real_type> &a){ v[0] -= a.v[0]; v[1] -= a.v[1]; return *this; }
	
	struct LexicalLess{ bool operator()(const TVec2<RealType> &v1, const TVec2<RealType> &v2) const{
		if(v1.v[0] < v2.v[0]){ return true; }
		else if(v1.v[0] == v2.v[0]){
			return v1.v[1] < v2.v[1];
		}else{
			return false;
		}
	}};
};

template <class RealType>
TVec2<RealType> operator * (const TVec2<RealType> &a, const RealType &c){
	TVec2<RealType> ret(a);
	ret *= c;
	return ret;
}
template <class RealType>
TVec2<RealType> operator * (const RealType &c, const TVec2<RealType> &a){
	TVec2<RealType> ret(a);
	ret *= c;
	return ret;
}
template <class RealType>
TVec2<RealType> operator / (const TVec2<RealType> &a, const RealType &c){
	TVec2<RealType> ret(a);
	ret /= c;
	return ret;
}
template <class RealType>
TVec2<RealType> operator + (const TVec2<RealType> &a, const TVec2<RealType> &b){
	TVec2<RealType> ret(a);
	ret += b;
	return ret;
}
template <class RealType>
TVec2<RealType> operator - (const TVec2<RealType> &a, const TVec2<RealType> &b){
	TVec2<RealType> ret(a);
	ret -= b;
	return ret;
}
template <class RealType>
TVec2<RealType> operator - (const TVec2<RealType> &a){
	return TVec2<RealType>(-a.v[0], -a.v[1]);
}

template <class RealType>
bool operator==(const TVec2<RealType> &A, const TVec2<RealType> &B){
	return ((A[0] == B[0]) && (A[1] == B[1]));
}
template <class RealType>
bool operator!=(const TVec2<RealType> &A, const TVec2<RealType> &B){
	return ((A[0] != B[0]) || (A[1] != B[1]));
}

template <class RealType>
std::ostream& operator<<(std::ostream &os, TVec2<RealType> &v){
	os << '{' << v.v[0] << ", " << v.v[1] << '}';
}

#endif // _TVEC2_H_
