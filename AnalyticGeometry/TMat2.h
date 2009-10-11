#ifndef _TMAT2_H_
#define _TMAT2_H_

#include <iostream>

template <class RealType>
class TMat2{
	typedef RealType real_type;
	real_type m[4]; // Column major
	enum special_ctor{
		CTOR_ZERO,
		CTOR_ID
	};
	TMat2(special_ctor s){
		switch(s){
		case CTOR_ZERO:
			for(size_t i = 0; i < 4; ++i){ m[i] = real_type(0); }
			break;
		case CTOR_ID:
			m[0] = m[3] = real_type(1);
			m[1] = m[2] = real_type(0);
			break;
		}
	}
public:
	static const TMat2<real_type> Zero;
	static const TMat2<real_type> Identity;
	
	TMat2(){
		for(size_t i = 0; i < 4; ++i){ m[i] = real_type(0); }
	}
	TMat2(const TMat2 &M){
		for(size_t i = 0; i < 4; ++i){ m[i] = M.m[i]; }
	}
	
	TMat2<real_type>& operator = (const TMat2<real_type> &mat){
		for(size_t i = 0; i < 4; ++i){ m[i] = mat.m[i]; } return *this;
	}

	TMat2<real_type>& operator+=(const TMat2<real_type>& mat){
		for(size_t i = 0; i < 4; ++i){ m[i] += mat.m[i]; } return *this;
	}
	TMat2<real_type>& operator-=(const TMat2<real_type>& mat){
		for(size_t i = 0; i < 4; ++i){ m[i] -= mat.m[i]; } return *this;
	}
	TMat2<real_type>& operator*=(const real_type &s){
		for(size_t i = 0; i < 4; ++i){ m[i] *= s; } return *this;
	}
	TMat2<real_type>& operator*=(const TMat2<real_type> &M){
		real_type n[2];
		n[0] = m[0]*M(0,0)+m[2]*M(1,0);
		n[1] = m[1]*M(0,0)+m[3]*M(1,0);
		m[2] = m[0]*M(0,1)+m[2]*M(1,1);
		m[3] = m[1]*M(0,1)+m[3]*M(1,1);
		m[0] = n[0];
		m[1] = n[1];
		return *this;
	}
	TMat2<real_type>& operator/=(const real_type &s){
		for(size_t i = 0; i < 4; ++i){ m[i] /= s; } return *this;
	}
	
	real_type& operator()(size_t row, size_t col){
		return m[2*col+row];
	}
	const real_type& operator()(size_t row, size_t col) const{
		return m[2*col+row];
	}

	TVec2<real_type> LinearTransform(const TVec2<real_type> &v) const{
		return TVec2<real_type>(v(0)*m[0] + v(1)*m[2], v(0)*m[1] + v(1)*m[3]);
	}
	real_type Trace() const{
		return m[0]+m[3];
	}
	real_type Det() const{
		return m[0]*m[3]-m[1]*m[2];
	}
	// Solve returns the determinant; can check for singularity
	real_type Solve(const TVec2<real_type> &RHS, TVec2<real_type> &Unknown) const{
		real_type det = Det();
		Unknown[0] = (m[3]*RHS[0]-m[2]*RHS[1])/det;
		Unknown[1] = (m[0]*RHS[1]-m[1]*RHS[0])/det;
		return det;
	}
	static TMat2<real_type> OuterProduct(const TVec2<real_type> &u, const TVec2<real_type> &v){
		TMat2<real_type> ret;
		for(size_t j = 0; j < 2; ++j){
			for(size_t i = 0; i < 2; ++i){
				ret(i,j) = u[i]*v[j];
			}
		}
		return ret;
	}
};

template <class RealType>
TVec2<RealType> operator*(const TMat2<RealType> &M, const TVec2<RealType> &v){
	return TVec2<RealType>(M(0,0)*v[0] + M(0,1)*v(1), M(1,0)*v(0) + M(1,1)*v(1));
}
template <class RealType>
TMat2<RealType> operator*(const RealType &scale, const TMat2<RealType> &M){
	TMat2<RealType> ret(M); ret *= scale;
	return ret;
}
template <class RealType>
TMat2<RealType> operator+(const TMat2<RealType> &A, const TMat2<RealType> &B){
	TMat2<RealType> ret(A); ret += B;
	return ret;
}
template <class RealType>
TMat2<RealType> operator-(const TMat2<RealType> &A, const TMat2<RealType> &B){
	TMat2<RealType> ret(A); ret -= B;
	return ret;
}
template <class RealType>
TMat2<RealType> operator-(const TMat2<RealType> &A){
	TMat2<RealType> ret(A);
	for(size_t j = 0; j < 2; ++j){
		for(size_t i = 0; i < 2; ++i){
			ret(i,j) = -ret(i,j);
		}
	}
	return ret;
}
template <class RealType>
TMat2<RealType> operator*(const TMat2<RealType> &A, const TMat2<RealType> &B){
	TMat2<RealType> ret(A); ret *= B;
	return ret;
}
template <class RealType>
bool operator==(const TMat2<RealType> &A, const TMat2<RealType> &B){
	return ((A(0,0) == B(0,0)) && (A(1,0) == B(1,0)) && (A(0,1) == B(0,1)) && (A(1,1) == B(1,1)));
}
template <class RealType>
bool operator!=(const TMat2<RealType> &A, const TMat2<RealType> &B){
	return ((A(0,0) != B(0,0)) || (A(1,0) != B(1,0)) || (A(0,1) != B(0,1)) || (A(1,1) == B(1,1)));
}

template <class RealType>
const TMat2<RealType> TMat2<RealType>::Zero(TMat2<RealType>::CTOR_ZERO);
template <class RealType>
const TMat2<RealType> TMat2<RealType>::Identity(TMat2<RealType>::CTOR_ID);

template <class RealType>
std::ostream& operator<<(std::ostream &os, const TMat2<RealType> &M){
	os << "{{" << M(0,0) << ", " << M(0,1) << "},{" << M(1,0) << ", " << M(1,1) << "}}";
}

#endif
