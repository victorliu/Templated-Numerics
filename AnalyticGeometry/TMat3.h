#ifndef _TMAT3_H_
#define _TMAT3_H_

template <class RealType>
class TMat3{
	typedef RealType real_type;
	real_type m[9]; // Column major
	
	static const TMat3 Identity;
	
	TMat3(){
		for(size_t i = 0; i < 9; ++i){ m[i] = 0; }
	}
	TMat3(const TMat3 &m){
		for(size_t i = 0; i < 9; ++i){ m[i] = m.m[i]; }
	}
	
	TMat3<real_type>& operator = (const TMat3<real_type> &mat){
		for(size_t i = 0; i < 9; ++i){ m[i] = mat.m[i]; } return *this;
	}

	TMat3<real_type>& operator+=(const TMat3<real_type>& mat){
		for(size_t i = 0; i < 9; ++i){ m[i] += mat.m[i]; } return *this;
	}
	TMat3<real_type>& operator-=(const TMat3<real_type>& mat){
		for(size_t i = 0; i < 9; ++i){ m[i] -= mat.m[i]; } return *this;
	}
	TMat3<real_type>& operator*=(const real_type &s){
		for(size_t i = 0; i < 9; ++i){ m[i] *= s; } return *this;
	}
	TMat3<real_type>& operator/=(const real_type &s){
		for(size_t i = 0; i < 9; ++i){ m[i] /= s; } return *this;
	}
	
	real_type& operator()(size_t row, size_t col){
		return m[3*col+row];
	}
	const real_type& operator()(size_t row, size_t col) const{
		return m[3*col+row];
	}
	
	TPt2<real_type> AffineTransform(const TPt2<real_type> &p) const{
		real_type ret[3];
		ret[0] = p(0)*m[0] + p(1)*m[3] + m[6];
		ret[1] = p(0)*m[1] + p(1)*m[4] + m[7];
		ret[2] = p(0)*m[2] + p(1)*m[5] + m[8];
		if(real_type::Zero == ret[2]){
			return TPt2<real_type>(ret[0], ret[1]);
		}else if(ret[2] < real_type::Zero){
			return TPt2<real_type>(-ret[0], -ret[1]);
		}else{
			return TPt2<real_type>(ret[0], ret[1]);
		}
	}
	TVec2<real_type> LinearTransform(const TVec2<real_type> &v) const{
		return TVec2<real_type>(v(0)*m[0] + v(1)*m[3], v(0)*m[1] + v(1)*m[4]);
	}
};

#endif
