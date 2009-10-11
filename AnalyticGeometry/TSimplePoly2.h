#ifndef _TSIMPLEPOLY2_H_
#define _TSIMPLEPOLY2_H_

#include <TVec2.h>
#include <TPt2.h>
#include <TMat3.h>

template <class RealType>
class TSimplePoly2{
public:
	typedef size_t size_type;
	typedef RealType real_type;
	typedef TPt2<real_type> point_type;
	typedef std::vector<point_type> pt_vec;
private:
	pt_vec pt;
	point_type min, max; // bounding box
	bool convex;
	
	void ComputeBoundingBox(){
		if(pt.size() < 1){ return; }
		min = pt[0];
		max = pt[0];
		for(pt_vec::const_iterator i = pt.begin();
		    pt.end() != i; ++i){
			if(min.r[0] > i->r[0]){ min.r[0] = i->r[0]; }
			if(min.r[1] > i->r[1]){ min.r[1] = i->r[1]; }
			if(max.r[0] < i->r[0]){ max.r[0] = i->r[0]; }
			if(max.r[1] < i->r[1]){ max.r[1] = i->r[1]; }
		}
		// Determine convexity
	}
public:
	TSimplePoly2(){}
	TSimplePoly2(const TSimplePoly2 &p):pt(p.pt),min(p.min),max(p.max),convex(p.convex){}
	TSimplePoly2(const pt_vec &pv):pt(pv){
		ComputeBoundingBox();
	}
	
	TSimplePoly2<real_type>& operator=(const TSimplePoly2<real_type> &p){
		pt = p.pt;
		min = p.min;
		max = p.max;
		convex = p.convex;
	}
	
	// If result is non simple, produces SOMETHING
	void Union(const TSimplePoly2<real_type> &other, TSimplePoly2<real_type> &result){
	}
	void Intersect(const TSimplePoly2<real_type> &other, TSimplePoly2<real_type> &result){
	}
	void Difference(const TSimplePoly2<real_type> &other, TSimplePoly2<real_type> &result){
	}
	void Xor(const TSimplePoly2<real_type> &other, TSimplePoly2<real_type> &result){
	}
	void Grow(const real_type &amount){
	}
	void Translate(const TVec2<real_type> &v){
		for(pt::const_iterator i = pt.begin(); i != pt.end(); ++i){
			*i += v;
		}
	}
	void Rotate(const real_type &degrees){
	}
	void Scale(const real_type &scale, const point_type &about = point_type::Origin){
	}
	void Transform(const TMat3<real_type> &m){
	}
	
	static void MakeTriangle(TSimplePoly2<real_type> &poly){}
	static void MakeSquare(TSimplePoly2<real_type> &poly){}
	static void MakeCircle(TSimplePoly2<real_type> &poly){}
	
	struct Triangulation{
		struct Triangle{
			size_t v[3];
		};
		std::vector<Triangle> triangles;
	};
	void Triangulate(Triangulation &triangulation) const{
	}
};

#endif // _TSIMPLEPOLY2_H_
