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
	typedef TVec2<real_type> Vec2;
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
	real_type Area() const{
		int n = pt.size();
		real_Type area(0);

		for(int p=n-1, q=0; q < n; p = q++){
			A += pt[p].r[0]*pt[q].r[1] - pt[q].r[0]*pt[p].r[1];
		}
		return (real_type(1)/real_type(2))*A;
	}
	
	static void MakeTriangle(TSimplePoly2<real_type> &poly){}
	static void MakeSquare(TSimplePoly2<real_type> &poly){}
	static void MakeCircle(TSimplePoly2<real_type> &poly){}
	
	struct Triangulation{
		struct Triangle{
			size_t v[3];
			Triangle(size_t a, size_t b, size_t c){
				v[0] = a; v[1] = b; v[2] = c;
			}
		};
		std::vector<Triangle> triangles;
	};
	void Triangulate(Triangulation &triangulation) const{
		const size_t n = pt.size();
		if(n < 3){ return; }
		// Make a copy of all the vertices
		std::vector<size_t> V; V.reserve(n);
		if(Area < real_type(0)){
			for(int v = 0; v < n; ++v){ V[v] = v; }
		}else{
			for(int v = 0; v < n; ++v){ V[v] = n-1-v; }
		}
		
		size_t nv = n;
		size_t count = 2*nv;
		for(size_t v = nv-1; nv > 2; ){
			if(0 >= (count--)){ return; } // bad polygon
			// get 3 consecutive vertices
			size_t u = v;   if(nv <= u){ u = 0; } // prev
			v = u+1;     if(nv <= v){ v = 0; } // mid
			size_t w = v+1; if(nv <= w){ w = 0; } // next
			
			// Can clip the ear?
			bool can_clip = true;
			do{
				Vec2 a(pt[V[v]] - pt[V[u]]);
				Vec2 b(pt[V[w]] - pt[V[u]]);
				if(Vec2::Cross(a,b) < 0){ can_clip = false; break; }
				for(size_t p = 0; p < nv; ++p){
					if((p == u) || (p == v) || (p == w)){ continue; }
					if(TTriangle<real_type>(pt[V[u]], a, b).Contains(pt[V[p]])){ can_clip = false; break; }
				}
			}while(false);
			
			// Clip off the ear
			if(can_clip){
				triangulation.triangles.push_back(V[u], V[v], V[w]);
				for(size_t s = v, t = v+1; t < nv; ++s, ++t){
					V[s] = V[t];
				}
				--nv;
				count = 2*nv;
			}
		}
	}
};

#endif // _TSIMPLEPOLY2_H_
