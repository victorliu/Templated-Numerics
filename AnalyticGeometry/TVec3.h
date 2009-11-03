#ifndef _TVEC3_H_
#define _TVEC3_H_

#include <iostream>

template <class RealType>
struct TVec3{
	typedef RealType value_type;
	typedef RealType real_type;
	real_type v[3];
	
	static const TVec3 Zero;
	
	TVec3(){ v[0] = 0; v[1] = 0; v[2] = 0; }
	TVec3(const real_type &x, const real_type &y, const real_type &z){
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}
	TVec3(const TVec3 &p){ v[0] = p.v[0]; v[1] = p.v[1]; v[2] = p.v[2]; }
	TVec3& operator=(const TVec3 &p){ v[0] = p.v[0]; v[1] = p.v[1]; v[2] = p.v[2]; return *this; }
	
	const real_type& operator[](size_t idx) const{ return v[idx]; }
	real_type& operator[](size_t idx){ return v[idx]; }
	const real_type& x() const{ return v[0]; }
	const real_type& y() const{ return v[1]; }
	const real_type& z() const{ return v[2]; }
	real_type& x(){ return v[0]; }
	real_type& y(){ return v[1]; }
	real_type& z(){ return v[2]; }

	void Normalize(){
		real_type L = Length();
		v[0] /= L;
		v[1] /= L;
		v[2] /= L;
	}

	/// Returns the L^2 norm of the vector.
	real_type Length() const{
		if(abs(v[1]) > abs(v[2])){
			if(abs(v[1]) > abs(v[0])){
				return v[1]*sqrt(real_type(1) + v[0]/v[1] + v[2]/v[1]);
			}else{
				return v[0]*sqrt(real_type(1) + v[1]/v[0] + v[2]/v[0]);
			}
		}else{
			if(abs(v[2]) > abs(v[0])){
				return v[2]*sqrt(real_type(1) + v[0]/v[2] + v[1]/v[2]);
			}else{
				return v[0]*sqrt(real_type(1) + v[1]/v[0] + v[2]/v[0]);
			}
		}
	}
	/// Returns the square of the L^2 norm of the vector.
	real_type LengthSq() const{
		return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	}

	/// Performs the dot (inner) product on two vectors (Euclidean R^2 metric).
	static real_type Dot(const TVec3<real_type> &a, const TVec3<real_type> &b){
		return a.v[0]*b.v[0] + a.v[1]*b.v[1] + a.v[2]*b.v[2];
	}
	
	/// Performs the scalar product (length of cross product).
	static TVec3<real_type> Cross(const TVec3<real_type> &a, const TVec3<real_type> &b){
		return TVec3<real_type>(
			a.v[1]*b.v[2] - a.v[2]*b.v[1],
			a.v[2]*b.v[0] - a.v[0]*b.v[2],
			a.v[0]*b.v[1] - a.v[1]*b.v[0]
			);
	}
	static TVec3<real_type> MakeTriad(const TVec3<real_type> &v0, TVec3<real_type> &v1, TVec3<real_type> &v2){
		TVec3<real_type> Given = v0;
		Given.Normalize();
		if(abs(Given.x) > abs(Given.y)){
			real_type invLen = real_type(1) / sqrt(Given.x*Given.x + Given.z*Given.z);
			v1 = TVec3<real_type>(-Given.z * invLen, 0, Given.x * invLen);
		}else{
			real_rype invLen = real_type(1) / sqrt(Given.y*Given.y + Given.z*Given.z);
			v1 = TVec3<real_type>(0, Given.z * invLen, -Given.y * invLen);
		}
		v2 = TVec3<real_type>::Cross(Given, v1);
	}

	TVec3<real_type>& operator *= (const real_type &C){ v[0] *= C; v[1] *= C; v[2] *= C; return *this; }
	TVec3<real_type>& operator /= (const real_type &C){ v[0] /= C; v[1] /= C; v[2] /= C; return *this; }
	TVec3<real_type>& operator += (const TVec3<real_type> &a){ v[0] += a.v[0]; v[1] += a.v[1]; v[2] += a.v[2]; return *this; }
	TVec3<real_type>& operator -= (const TVec3<real_type> &a){ v[0] -= a.v[0]; v[1] -= a.v[1]; v[2] -= a.v[2]; return *this; }
};

template <class RealType>
TVec3<RealType> operator * (const TVec3<RealType> &a, const RealType &c){
	TVec3<RealType> ret(a);
	ret *= c;
	return ret;
}
template <class RealType>
TVec3<RealType> operator * (const RealType &c, const TVec3<RealType> &a){
	TVec3<RealType> ret(a);
	ret *= c;
	return ret;
}
template <class RealType>
TVec3<RealType> operator / (const TVec3<RealType> &a, const RealType &c){
	TVec3<RealType> ret(a);
	ret /= c;
	return ret;
}
template <class RealType>
TVec3<RealType> operator + (const TVec3<RealType> &a, const TVec3<RealType> &b){
	TVec3<RealType> ret(a);
	ret += b;
	return ret;
}
template <class RealType>
TVec3<RealType> operator - (const TVec3<RealType> &a, const TVec3<RealType> &b){
	TVec3<RealType> ret(a);
	ret -= b;
	return ret;
}
template <class RealType>
TVec3<RealType> operator - (const TVec2<RealType> &a){
	return TVec3<RealType>(-a.v[0], -a.v[1]);
}

template <class RealType>
bool operator==(const TVec3<RealType> &A, const TVec3<RealType> &B){
	return ((A[0] == B[0]) && (A[1] == B[1]) && (A[2] == B[2]));
}
template <class RealType>
bool operator!=(const TVec3<RealType> &A, const TVec3<RealType> &B){
	return ((A[0] != B[0]) || (A[1] != B[1]) || (A[2] != B[2]));
}

template <class RealType>
std::ostream& operator<<(std::ostream &os, TVec3<RealType> &v){
	os << '{' << v.v[0] << ", " << v.v[1] << ", " << v.v[2] << '}';
}



template <class RealType>
struct TVec3LexicalLess{ bool operator()(const TVec3<RealType> &v1, const TVec3<RealType> &v2) const{
	if(v1.v[0] < v2.v[0]){ return true; }
	else if(v1.v[0] == v2.v[0]){
		if(v1.v[1] < v2.v[1]){ return true; }
		else if(v1.v[1] == v2.v[1]){ return v1.v[2] < v2.v[2]; }
		else{ return false; }
	}else{ return false; }
}};

#endif // _TVEC3_H_
