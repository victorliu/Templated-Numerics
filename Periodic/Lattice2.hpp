#ifndef _LATTICE_HPP_
#define _LATTICE_HPP_

#include <AnalyticGeometry/TVec2.h>
#include <AnalyticGeometry/TPt2.h>
#include <AnalyticGeometry/TMat2.h>
#include <vector>
#include <map>

/* Lattice2 represents a 2D point lattice.
 * The real number type must must be a subfield of the real number field
 * and support exact equality comparison (a == b must always be true), which
 * rules out floating point types and fixed point types (due to round off
 * error). Good choices are rationals and rationals adjoined with radicals.
 * Efforts were made to ensure all operations needed are closed in those
 * fields.
 */

template <class RealType>
class Lattice2{
public:
	typedef RealType real_type;
	typedef TVec2<real_type> Vec2;
	typedef TPt2<real_type> Pt2;
	typedef TVec2<int> Coord2;
	typedef TVec2<real_type> UVCoord;
	typedef TMat2<real_type> Mat2;
	typedef TVec2<real_type> RecipVec2;
	// In 2D, there are only 5 distinct lattice symmetry groups
	enum LatticeSymmetry{
		OBLIQUE,
		RHOMBIC,
		RECTANGULAR,
		SQUARE,
		HEXAGONAL,
		DEGENERATE
	};
	
	// We provide no default constructor since that would be meaningless.
	// We provide no copy constructor since these objects are so light.
	
	// It is assumed that basis1 and basis2 are not collinear (that
	// they span the space). No explicit check of this is ever
	// performed.
	Lattice2(const Vec2 &basis1, const Vec2 &basis2):u(basis1),v(basis2){
		Reduce(u,v);
		if(Vec2::Cross(u,v) < 0){
			std::swap(u,v);
		}
	}
	
	// Returns the (reduced) lattice basis. Note that this may not be the
	// same as the basis provided to the constructor even though it was
	// already reduced.
	void GetBasis(Vec2 &basis1, Vec2 &basis2) const{ basis1 = u; basis2 = v; }
	const Vec2& Basis(size_t which) const{ return (0 == which) ? u : v; }
	
	 // Scaled by down 2pi, as in, the actual reciprocal basis is
	 // 2*M_PI*g1, 2*M_PI*g2.
	void GetReciprocalBasis(RecipVec2 &g1, RecipVec2 &g2) const{
		real_type uv = Vec2::Cross(u,v);
		g1 = Vec2::Rotate90CCW(v) / (-uv);
		g2 = Vec2::Rotate90CCW(u) /   uv;
	}
	void Reciprocate(){ // flip primal and reciprocal basis
		Vec2 up, vp;
		GetReciprocalBasis(up, vp);
		u = up;
		v = vp;
	}
	
	// This is a function specially tailored for 2D symmetries, but it
	// should still generalize well. Part of me thinks it is inelegant
	// to special case the symmetry "detection" for each lattice type,
	// but a general symmetry search is quite complicated.
	LatticeSymmetry GetSymmetries(std::vector<Mat2> &symmetries) const{
		/*
		Oblique:      (3,0) (1,2)
		Square:       (1,0) (0,1)               u.v = 0, |u| = |v|
		Rectangular:  (1,0) (0,2)               u.v = 0
		Hexagonal:    (1,0) (1/2,sqrt(3)/2)    |u.v| = 1/2, |u| = |v|
		Rhombic:      (1,0) (1/2, 1)           |u.v| = 1/2 or equivalently |u| = |v|

				  |             |u.v|
				  +----------------------------------
		|u|==|v|  |   else   |    0    |     1/2     
		----------+----------+---------+------------
		   no     | oblique  | rect    |  rhombic
		   yes    | rhombic  | square  |  hex        
		*/
		
		if(real_type(0) == Vec2::Cross(u,v)){
			return DEGENERATE;
		}
		
		symmetries.clear();
		symmetries.push_back(Mat2::Identity);  // These two are always symmetries
		symmetries.push_back(-Mat2::Identity);
		real_type u2 = u.LengthSq();
		real_type v2 = v.LengthSq();
		real_type uv = abs(Vec2::Dot(u,v));
		real_type hu = real_type(2)*uv - u2; // 2*|u.v| - |u|^2, equal to zero if |u.v|/|u| == 1/2 |u|
		real_type hv = real_type(2)*uv - v2;
		
		int row = 0, col = 0; // index into the table above
		if(u2 == v2){ row = 1; }
		bool udiag = false;
		if(uv == 0){
			col = 1;
		}else if(hu == 0){
			udiag = true;
			col = 2;
		}else if(hv == 0){
			col = 2;
		}
		
		switch(col*2+row){
		case 0:
			return OBLIQUE;
		case 1: // two sides of rhombus
			symmetries.push_back(FlipAlong(u+v));
			symmetries.push_back(FlipAlong(u-v));
			return RHOMBIC;
		case 4: // side and a diagonal of rhombus
			if(udiag){
				symmetries.push_back(FlipAlong(u));
				symmetries.push_back(FlipAlong(Vec2::Rotate90CCW(u)));
			}else{
				symmetries.push_back(FlipAlong(v));
				symmetries.push_back(FlipAlong(Vec2::Rotate90CCW(v)));
			}
			return RHOMBIC;
		case 2:
			symmetries.push_back(FlipAlong(u));
			symmetries.push_back(FlipAlong(v));
			return RECTANGULAR;
		case 3:
			symmetries.push_back(Rot90n(1));
			//symmetries.push_back(Rot90n(2)); // done above
			symmetries.push_back(Rot90n(3));
			symmetries.push_back(FlipAlong(u));
			symmetries.push_back(FlipAlong(v));
			symmetries.push_back(FlipAlong(u+v));
			symmetries.push_back(FlipAlong(u-v));
			return SQUARE;
		case 5:
			symmetries.push_back(Rot60n(1));
			symmetries.push_back(Rot60n(2));
			//symmetries.push_back(Rot60n(3)); // done above
			symmetries.push_back(Rot60n(4));
			symmetries.push_back(Rot60n(5));
			symmetries.push_back(FlipAlong(u));
			symmetries.push_back(FlipAlong(v));
			symmetries.push_back(FlipAlong(u+v));
			symmetries.push_back(FlipAlong(u-v));
			symmetries.push_back(FlipAlong(Vec2::Rotate90CCW(u)));
			symmetries.push_back(FlipAlong(Vec2::Rotate90CCW(v)));
			return HEXAGONAL;
		}
	}
	
protected:
	Vec2 u, v;
	
	// Generators of rigid transformations
	static Mat2 FlipAlong(const Vec2 &u){
		return Mat2::Identity - (real_type(2)/u.LengthSq())*Mat2::OuterProduct(u,u);
	}
	static Mat2 Rot90n(int n){
		static const int costab[4] = {1, 0, -1, 0};
		static const int sintab[4] = {0, 1, 0, -1};
		Mat2 rot;
		rot(1,1) = rot(0,0) = costab[n%4];
		rot(0,1) = -(rot(1,0) = sintab[n%4]);
		return rot;
	}
	static Mat2 Rot60n(int n){
		static const real_type half(real_type(1)/real_type(2));
		static const real_type r3 = sqrt((real_type(3)/real_type(4)));
		static const real_type costab[6] = {1, half, -half, -1, -half, half};
		static const real_type sintab[6] = {0, r3, r3, 0, -r3, -r3};
		Mat2 rot;
		rot(1,1) = rot(0,0) = costab[n%6];
		rot(0,1) = -(rot(1,0) = sintab[n%6]);
		return rot;
	}
	
	// Basis reduction
	// Implements Lagrange's algorithm
	static void Reduce(Vec2 &u, Vec2 &v){
		if(u.LengthSq() > v.LengthSq()){
			std::swap(u,v);
		}

		Mat2 G; // Gram matrix; G(1,0) is never used
		G(0,0) = Vec2::Dot(u,u);
		G(0,1) = Vec2::Dot(u,v);
		G(1,1) = Vec2::Dot(v,v);
		
		do{
			real_type x = round(G(0,1) / G(0,0));
			real_type y = G(0,1) - x*G(0,0); // remainder of the above division
			Vec2 r = v - x*u;
			v = u;
			u = r;
			// Update Gram matrix
			std::swap(G(0,0), G(1,1));
			G(0,0) -= x*(y+G(0,1));
			G(0,1) = y;
		}while(u.LengthSq() < v.LengthSq());
	}
	
private:
	// computes floor(sqrt(u)) by Newton iteration
	static int isqrt(int u){
		if(u <= 0){ return 0; }
		int a = 2;
		int b = u/a;
		while(((a-b) > 1) || ((b-a) > 1)){
			a = (a+b)/2;
			b = u/a;
		}
		return (a<b)?a:b;
	}
	// computes ceil(sqrt(u))
	static int isqrt_ceil(int u){
		int s = isqrt(u);
		return (s*s == u) ? s : s+1;
	}
public:
	/* A ring is a set of vectors S_m = {x_1, ..., x_N}, where |R_m| = N
	 * such that |x_i| = R_m, and the rings are ordered such that R_i < R_{i+1}
	 * In other words, the rings comprise sets of points with equal magnitude.
	 * Note that points in a ring do not necessarily correspond to a single
	 * symmetry group of the lattice. A simple counterexample is for a square
	 * lattice where (5,5) and (1,7) have the same radius, but are not related
	 * any symmetry of the square lattice.
	 *
	 * The points vector will contain all the points on the ring.
	 * Runs in time O(R_m), i.e. in time proportional to the perimeter of the
	 * ring.
	 */
	void GetNextPointRing(int &a, int &b, std::vector<Coord2> &points) const{
		/*
		next = x = k*u + l*v, and need |x| > r = |current|
		So find smallest k such that
		  k^2 |u|^2 + l^2 |v|^2 + 2 l k v.u == r^2
		has no solution for any l. Occurs for
		               r^2 |v|^2
		  k^2 > ----------------------- = t
		         |u|^2 |v|^2 - (u.v)^2
		
		k = 1+floor(sqrt(t)) = 1+floor(sqrt(floor(t)))
		*/
		real_type u2 = u.LengthSq(), v2 = v.LengthSq(), uv = Vec2::Dot(u,v);
		real_type r2 = (real_type(a)*u + real_type(b)*v).LengthSq();
		real_type uv_v2 = uv/v2;
		real_type t = r2*v2 / (u2*v2 - uv*uv);
		
		int kmax = 1+isqrt((int)floor(t));
		real_type next_r2 = r2 + u2 + u2 + v2 + v2;
		
		for(int k = -kmax; k <= kmax; ++k){
			real_type k2(k*k);
			real_type alpha = real_type(-k)*uv_v2;
			real_type beta  = (k2*u2-r2)/v2;
			// roots are alpha +/- sqrt(alpha^2 - beta)
			int l[4]; // there are always at most 4 possible values of l to check
			int nl = 4;
			// Now figure out what nearest l's for current k
			// check discriminant to see if there is a solution
			real_type d = alpha*alpha-beta;
			if(d < 0){ // no solution, find the minimizer instead
				l[0] = (int)floor(alpha); // minimizer is alpha rounded to nearest int
				l[1] = l[0] + 1;     // so just check the floor and ceil, this may be excessive
				if(l[0] == 0){ // check for duplication
					l[2] = -l[1];
					nl = 3;
				}else if(l[1] == 0){
					l[2] = -l[0];
					nl = 3;
				}else{
					l[2] = -l[0];
					l[3] = -l[1];
				}
			}else if(d == 0){ // exactly on the circle
				l[0] = (int)floor(alpha);
				l[1] = l[0] + 1;
				if(l[0] == 0){ // check for duplication
					l[2] = -l[1];
					nl = 3;
				}else if(l[1] == 0){
					l[2] = -l[0];
					nl = 3;
				}else{
					l[2] = -l[0];
					l[3] = -l[1];
				}
			}else{
				// we really want to the roots rounded away from the average of the roots:
				//   floor(alpha - sqrt(d))   and   ceil(alpha + sqrt(d))
				// We note that
				//   floor(alpha - sqrt(d)) = floor(alpha) + floor(-sqrt(d)) + r, where r = 0 or 1
				//                          = floor(alpha) - ceil(sqrt(d)) + r
				//                          = floor(alpha) - isqrt_ceil(ceil(d)) + r
				// Similarly
				//   ceil(alpha + sqrt(d)) = ceil(alpha) + ceil(sqrt(d)) - r
				//                         = ceil(alpha) + isqrt_ceil(ceil(d)) - r
				int ic = isqrt_ceil((int)ceil(d));
				l[0] = (int)floor(alpha) - ic;
				l[1] = l[0] + 1;
				l[2] = (int)ceil(alpha) + ic;
				if(l[2] == l[0]){ // check for duplication
					l[2] = l[2] - 1;
					nl = 3;
				}else{
					l[3] = l[2] - 1;
				}
			}
			for(int i = 0; i < nl; ++i){
				Vec2 x = real_type(k)*u + real_type(l[i])*v;
				real_type x2 = x.LengthSq();
				if(x2 > r2 && x2 < next_r2){
					points.clear(); // if we found a smaller radius, get rid of all the old points
					next_r2 = x2;
					a = k; b = l[i];
					points.push_back(Coord2(k,l[i]));
				}else if(x2 == next_r2){
					points.push_back(Coord2(k,l[i]));
				}
			}
		}
	}
	
	// Here are the Voronoi region related functions.
	// If the lattice is a real space lattice, the Voronoi region is the
	// Wigner-Seitz cell. If the lattice is a reciprocal space lattice,
	// then the Voronoi region is the first Brillouin zone.
	
	// points will contain 3 points, except when the lattice is
	// square or rectangular, when it will contain 2. Only half
	// the points are included, since the negatives of the ones
	// returned are implicit.
	// This function assumes that u and v are the shortest possible
	// basis vectors (which is something that Reduce should guarantee).
	void GetVoronoiDefiningPoints(std::vector<Pt2> &points) const{
		points.clear();
		points.push_back(Pt2::Origin + u);
		//points.push_back(-u);
		points.push_back(Pt2::Origin + v);
		//points.push_back(-v);
		Vec2 upv(u+v), umv(u-v);
		real_type upv2 = upv.LengthSq(), umv2 = umv.LengthSq();
		if(umv2 < upv2){
			points.push_back(Pt2::Origin + umv);
			//points.push_back(-umv);
		}else if(upv2 < umv2){
			points.push_back(Pt2::Origin + upv);
			//points.push_back(-upv);
		}
	}
	// Given an arbitrary vector, fold it into the Voronoi region
	// by subtracting integral multiples of the basis vectors.
	void FoldIntoVoronoi(const std::vector<Pt2> &vpts, Pt2 &a) const{
		for(typename std::vector<Pt2>::const_iterator i = vpts.begin();i != vpts.end(); ++i){
			const Vec2 &vp = *i - Pt2::Origin; // current Voronoi defining point
			real_type vp2 = vp.LengthSq();
			// Each Voronoi defining point defines a hyperplane that is
			// equidistant between it and the origin.
			// We check if the vector is on the side closer to the origin.
			// If not, subtract a proper (integral) amount of vp.
			
			// We want to see if a projected onto vp is longer than vp/2.
			real_type d = Vec2::Dot(a-Pt2::Origin,vp)/vp2;
			real_type half(real_type(1)/real_type(2));
			int n = floor(d + half);
			a -= real_type(n)*vp;
		}
	}
	void FoldIntoVoronoi(Vec2 &a) const{
		std::vector<Pt2> vpts;
		GetVoronoiDefiningPoints(vpts);
		FoldIntoVoronoi(vpts, a);
	}
	Mat2 GramMatrix() const{
		Mat2 ret;
		ret(0,0) = u.LengthSq();
		ret(0,1) = ret(1,0) = Vec2::Dot(u,v);
		ret(1,1) = v.LengthSq();
		return ret;
	}
	Vec2 UVCoords(Pt2 &a) const{
		Mat2 T;
		T(0,0) = u.v[0];
		T(1,0) = u.v[1];
		T(0,1) = v.v[0];
		T(1,1) = v.v[1];
		Vec2 st; // coeffs of a in u,v basis
		T.Solve(a - Pt2::Origin, st);
		return st;
	}

#ifdef USING_NUMERIC_TYPE_TRAITS
	template <class OtherNumericType>
	void UVCoords(OtherNumericType &x, OtherNumericType &y) const{
		OtherNumericType xx(x), yy(y);
		OtherNumericType d(ScalarTraits<real_type>::numeric_value<OtherNumericType>(Vec2::Cross(u,v)));
		x = (ScalarTraits<real_type>::numeric_value<OtherNumericType>(v.v[1])*xx - ScalarTraits<real_type>::numeric_value<OtherNumericType>(v.v[0])*yy) / d;
		y = (ScalarTraits<real_type>::numeric_value<OtherNumericType>(u.v[0])*yy - ScalarTraits<real_type>::numeric_value<OtherNumericType>(u.v[1])*xx) / d;
	}
#endif

	// Convert from coefficients to a point (opposite of UVCoords)
	Pt2 operator()(const real_type &s, const real_type &t) const{
		return Pt2(Pt2::Origin + s*u + t*v);
	}
	Pt2 operator()(const UVCoord &st) const{
		return Pt2(Pt2::Origin + st[0]*u + st[1]);
	}
	Pt2 operator()(int s, int t) const{
		return Pt2(Pt2::Origin + real_type(s)*u + real_type(t)*v);
	}
	
	void FoldIntoParallelogram(Pt2 &a) const{
		Vec2 st(UVCoords(a));
		st[0] = st[0] - floor(st[0]);
		st[1] = st[1] - floor(st[1]);
		a = (*this)(st);
	}
	
	// Returns the volume of the fundamental parallelogram.
	// Also equal to the volume of the voronoi region.
	real_type Volume() const{
		return Vec2::Cross(u,v);
	}
	
	// Determine if a point is within the Voronoi region.
	bool InVoronoi(const std::vector<Pt2> &vpts, const Pt2 &a) const{
		for(typename std::vector<Pt2>::const_iterator i = vpts.begin();i != vpts.end(); ++i){
			const Vec2 &vp = *i - Pt2::Origin; // current Voronoi defining point
			real_type vp2 = vp.LengthSq();
			// Each Voronoi defining point defines a hyperplane that is
			// equidistant between it and the origin.
			// We check if the vector is on the side closer to the origin.
			
			// We want to see if a projected onto vp is longer than vp/2.
			real_type d = Vec2::Dot(real_type(2)*a,vp)/vp2;
			if(abs(d) > real_type(1)){ return false; }
		}
		return true;
	}
	bool InVoronoi(const Pt2 &a) const{
		std::vector<Pt2> vpts;
		GetVoronoiDefiningPoints(vpts);
		InVoronoi(vpts, a);
	}
	
	// The irreducible part of the Voronoi region is the part that can
	// completely regenerate all other points in the Voronoi region through
	// symmetry operations of the lattice. For example, for a square lattice,
	// it is the one-eigth of Voronoi region that is an isoceles triangle.
	// The irreducible part is determined by the rotational symmetries of
	// the lattice, and is always a wedge centered at the origin.
	//
	// Here, we return vectors which define normal vectors of halfspaces.
	// A point is in the irreducible part if the dot product with all
	// three vectors is positive. (Of course, must also be in the Voronoi
	// region, which is not part of this function.)
	// For less symmetric lattices, these vectors may not be unique.
	void GetIrreducibleVoronoiWedge(Vec2 &n1, Vec2 &n2, Vec2 &n3) const{
		// TODO
	}
	
	// Here begins code for "special points" methods of integration
	// over a unit cell of the lattice. The basic idea is to pick
	// points in the unit cell at which the function to be integrated
	// is to be evaluated, as well as the weighting factors to be
	// used to sum them up. These methods try to minimize the number
	// of points needed while cancelling off as many low order
	// Fourier components as possible.
public:
	struct SpecialPoint{
		Pt2 x;
		real_type weight;
		SpecialPoint(const Pt2 &v, const real_type &w):x(v),weight(w){}
	};
private:
public:
	/* Chadi-Cohen points are recursively determined within the first BZ. It
	 * can be thought of as a multi-dimensional generalization of the periodic
	 * midpoint rule.
	 */
	/*
	// It turns out that for the special point determination depends on the
	// symmetry of the lattice, and so should be special cased for each lattice
	void GetChadiCohenPoints(int iterations, std::vector<SpecialPoint> &points) const{
		// bootstrap; find the first point
		std::vector<Vec2> R_points;
		int cu,cv;
		GetNextPointRing(cu,cv, R_points);
		Vec2 k1 = R_points[0]/(4*R_points[0].LengthSq());
		
		for(int i = 0; i < iterations; ++i){
			// Find next points that are not zeroed out
			for(size_t ri = 0; ri < R_points.size(); ++ri){
				
			}
		}
	}
	*/
	/* Monkhorst-Pack points are uniformly distributed in the fundamental
	 * parallelogram of the lattice. They need to be folded back and reduced
	 * by the symmetry operations of the lattice. (We prefer points within
	 * the first BZ, which the original MP formulation does not care about.
	 *
	 * q is the sampling density in each dimension; there should be q^2 points
	 * counting multiplicities.
	 *
	 * Currently points related by symmetry operations of the lattice are not
	 * combined (doing so would result in points in only the irreducible BZ).
	 */
	void GetMonkhorstPackPoints(int q, std::vector<SpecialPoint> &points) const{
		points.clear();
		// First find the BZ
		std::vector<Pt2> vpts;
		GetVoronoiDefiningPoints(vpts);
		
		typedef std::map<Pt2,int,TPt2LexicalLess<typename Pt2::real_type> > vec_mult_map;
		vec_mult_map cand_pts;
		for(int i = 0; i < q; ++i){
			real_type ui(real_type(1+2*i-q)/real_type(2*q));
			for(int j = 0; j < q; ++j){
				real_type uj(real_type(1+2*j-q)/real_type(2*q));
				Pt2 x = Pt2::Origin + ui*u + uj*v;
				FoldIntoVoronoi(vpts, x);
				cand_pts[x]++;
			}
		}
		real_type iq2(real_type(1)/real_type(q*q));
		for(typename vec_mult_map::const_iterator it = cand_pts.begin(); it != cand_pts.end(); ++it){
			points.push_back(SpecialPoint(it->first, real_type(it->second)*iq2));
		}
	}
	// Same as above, but multiplies q by 3 (or an odd number) to generate a set of refined points
	void GetMonkhorstPackPointsRefinement(int &q, std::vector<SpecialPoint> &points, int refinement = 3) const{
		points.clear();
		// First find the BZ
		std::vector<Pt2> vpts;
		GetVoronoiDefiningPoints(vpts);
		
		if(0 == (refinement&1)){ ++refinement; }
		q *= refinement;
		
		typedef std::map<Pt2,int,TPt2LexicalLess<typename Pt2::real_type> > vec_mult_map;
		vec_mult_map cand_pts;
		for(int i = 0; i < q; ++i){
			// if refinement*(1+2*p) == 1+2*i for an integral p, then this point was already gotten
			//    refinement*2*p == 1+2*i-refinement
			if(0 == ((1+2*i-refinement)%(2*refinement))){ continue; }
			real_type ui(real_type(1+2*i-q)/real_type(2*q));
			for(int j = 0; j < q; ++j){
				if(0 == ((1+2*j-refinement)%(2*refinement))){ continue; }
				real_type uj(real_type(1+2*j-q)/real_type(2*q));
				Pt2 x = Pt2::Origin + ui*u + uj*v;
				FoldIntoVoronoi(vpts, x);
				cand_pts[x]++;
			}
		}
		real_type iq2(real_type(1)/real_type(q*q));
		for(typename vec_mult_map::const_iterator it = cand_pts.begin(); it != cand_pts.end(); ++it){
			points.push_back(SpecialPoint(it->first, real_type(it->second)*iq2));
		}
	}
	
	
	template <class NumericDomainType, class NumericFunctionType>
	class Function{
		virtual NumericFunctionType operator()(const TPt2<NumericDomainType> &x) const = 0;
		NumericFunctionType EvalVoronoi(TPt2<NumericDomainType> x, const std::vector<Pt2> vpts) const{
			FoldIntoVoronoi(vpts, x);
			return (*this)(x);
		}
	};
	template <class NumericType>
	struct IntegrationParameters{
		enum Method{
			MONKHORST_PACK,
			ADAPTIVE_CUBATURE
		};
		Method method;
		int max_num_points;
		NumericType abs_err, rel_err;
	};
	template <class NumericDomainType, class NumericFunctionType>
	NumericFunctionType Integrate(const Function<NumericDomainType,NumericFunctionType> &func, const IntegrationParameters<NumericFunctionType> &params) const{
		int q = isqrt(params.max_num_points);
		if(0 == (q&1)){ --q; }
		if(0 == q){ q = 1; }
		
		std::vector<SpecialPoint> points;
		GetMonkhorstPackPoints(q, points);
		
		//
	}
};

#endif
