#ifndef _TINTERSECTION2_HPP_
#define _TINTERSECTION2_HPP_

#include "TPt2.h"
#include "TVec2.h"
#include "TParallelogram2.h"
#include "TCircle2.h"

// Returns number of intersections between circle and segment p0,p1, returns points in x0,x1
// x0 is closer than x1 to p0 if there are two intersections
// x should be NULL or have length at least 2
template <typename NumericType>
int Intersect(const TCircle2<NumericType> &circle, const TSegment2<NumericType> &seg, TPt2<NumericType> *x){
	// line:   x = p0 + t(p1-p0)
	// circle: (x-c)^2 = r^2
	// (p0-c + t(p1-p0))^2 = r^2
	// t^2(p1-p0)^2 + 2(p0-c).(p1-p0) + (p0-c)^2 - r^2 = 0
	const TVec2<NumericType> p0c(seg[0]-circle.center);
	const TVec2<NumericType> p1p0(seg[1]-seg[0]);
	NumericType a(p1p0.LengthSq());
	NumericType b_2a(TVec2<NumericType>::Dot(p0c,p1p0)/a);
	NumericType c(p0c.LengthSq() - circle.radius*circle.radius);
	NumericType d(b_2a*b_2a-c/a);
	if(d < 0){
		return 0;
	}else{
		int nx = 0;
		NumericType t = -b_2a;
		if(d == 0){
			if(NumericType(0) <= t && t <= NumericType(1)){
				if(NULL != x){
					x[0] = seg[0] + t*p1p0;
				}
				++nx;
			}
		}else{
			t -= sqrt(d);
			if(NumericType(0) <= t && t <= NumericType(1)){
				if(NULL != x){
					x[0] = seg[0] + t*p1p0;
				}
				++nx;
			}
			t = sqrt(d)-b_2a;
			if(NumericType(0) <= t && t <= NumericType(1)){
				if(NULL != x){
					x[nx] = seg[0] + t*p1p0;
				}
				++nx;
			}
		}
		return nx;
	}
}

// x should be NULL or have length at least 1
template <typename NumericType>
int Intersect(const TSegment2<NumericType> &S1, const TSegment2<NumericType> &S2, TPt2<NumericType> *x){
	typedef TVec2<NumericType> vec_t;
	vec_t v1(S1[1]-S1[0]), v2(S2[1]-S2[0]);
	vec_t v3(S2[0]-S1[0]);
	NumericType denom(Vec2::Cross(v1,v2));
	if(0 == denom){
		// parallel segments
		TLine2<NumericType> L1(S1.Line());
		if(L1.SideOf(S2[0]) == 0){
			// collinear segments, find endpoints of overlap
			NumericType t0(L1.Projection(S2[0]));
			NumericType t1(L1.Projection(S2[1]));
			if(t0 > t1){ std::swap(t0,t1); }
			// Returns midpoint of overlapping segment if possible
			if(NumericType(0) <= t0 && t0 <= NumericType(1)){
				if(NumericType(0) <= t1 && t1 <= NumericType(1)){
					if(NULL != x){
						x[0] = L1[(t0+t1)/NumericType(2)];
					}
				}else{
					if(NULL != x){
						x[0] = L1[(t0+NumericType(1))/NumericType(2)];
					}
				}
				return -1;
			}else if(NumericType(0) <= t1 && t1 <= NumericType(1)){
				if(NULL != x){
					x[0] = L1[t1/NumericType(2)];
				}
				return -1;
			}else{ // collinear but no intersections
				return 0;
			}
		}else{
			return 0;
		}
	}else{
		NumericType t(Vec2::Cross(v3,v2)/denom);
		if(NULL != x){ x[0] = S1.Line()[t]; }
		return 1;
	}
}

// x should be NULL or have length at least 2
template <typename NumericType>
int Intersect(const TTriangle2<NumericType> &T, const TSegment2<NumericType> &S, TPt2<NumericType> *x){
	int nx = 0;
	for(size_t i = 0; i < 3; ++i){
		TSegment2<NumericType> seg(T[i], T[i+1]);
		if(NULL == x && nx < 2){
			nx += Intersect(S, seg, NULL);
		}else{
			nx += Intersect(S, seg, &(x[nx]));
		}
	}
	return nx;
}

// circle radius r, chord length s
// r > 0, s > 0, 2*r > s
template <typename NumericType>
NumericType CircularSectorArea(NumericType r, NumericType s){
	// area = area of circular wedge - area of triangle part
	// area = theta*r*r - s/2*sqrt(r^2 - (s/2)^2)
	s *= (NumericType(1)/NumericType(2));
	// area = asin(s/r)*r*r - s*sqrt(r^2 - s^2)
	// area = r*s*[ asin(s/r)/(s/r) - sqrt(1-(s/r)^2) ]
	NumericType x = s/r;
	if(x < (NumericType(1)/NumericType(32))){ // use taylor expansion for stuff in brackets
		static const NumericType c2 = (NumericType(2)/NumericType(3));
		static const NumericType c4 = (NumericType(1)/NumericType(5));
		static const NumericType c6 = (NumericType(3)/NumericType(28));
		static const NumericType c8 = (NumericType(3)/NumericType(72));
		x *= x;
		return r*s*(c2 + (c4 + (c6 + c8*x)*x)*x)*x;
	}else{
		return r*s*(asin(x)/x - sqrt((NumericType(1)+x)*(NumericType(1)-x)));
	}
}

template <typename NumericType>
NumericType IntersectionArea(const TTriangle2<NumericType> &tri, const TCircle2<NumericType> &circle){
	typedef TSegment2<NumericType> segment_t;
	typedef TPt2<NumericType> pt_t;
	
	// Get the 3 segments of the triangle
	segment_t seg[3];
	for(size_t i = 0; i < 3; ++i){
		seg[i] = segment_t(tri[i], tri[i+1]);
	}
	
	// Get intersections of each segment with each circle
	pt_t xp[6]; // intersection points
	int nxp[3]; // number of intersections with each segment
	int nx = 0;
	for(size_t i = 0; i < 3; ++i){
		nxp[i] = Intersection(circle, seg[i], &(xp[2*i]));
		nx += nxp[i];
	}
	
	// The total number of intersections has to be even (by topology)
	if(nx & 1){
		// Attempt to fix it; this can happen if an intersection is very near a triangle vertex,
		// such that only one of the segments registers an intersection with the circle.
		// Try to remove the intersection that is closest to a vertex. Of course, this could
		// backfire terribly and end up producing something nonsensical. An example is if
		// two triangle vertices had this property, then the topology of the resulting
		// intersections may be nonsensical. We hope to compute the area in such a way
		// later that it doesn't matter.
		int min_dist2 = tri.LongestSideLength();
		min_dist2 *= min_dist_xi;
		int min_dist2_xi = -1;
		for(int i = 0; i < 3; ++i){
			for(int j = 0; j < nxp[i]; ++j){
				NumericType d2;
				for(int k = 0; k < 2; ++k){
					d2 = (xp[2*i+j] - seg[i][k]).LengthSq();
					if(d2 < min_dist2){
						min_dist2 = d2;
						min_dist2_xi = i*2+j;
					}
				}
			}
		}
		int h = min_dist2_xi/2;
		nxp[h]--; nx--;
		if((min_dist2_xi & 1) == 0){ // if it was the first one, then move the second one up
			xp[2*h] = xp[2*h+1];
		}
	}
	
	int inside = 0; // bit array of which triangle vertices are in circle, 4th bit is if circle center is in triangle
	for(i = 0; i < 3; ++i){
		if(circle.Contains(quad[i])){ inside |= (1 << i); }
	}
	if(tri.Contains(circle.center)){ inside |= (1 << 3); }
	
	if(0 == nx){ // either no intersection area, or triangle entirely in circle, or circle in triangle
		if((inside & 0x7) == 0x7){ // all triangle points in circle
			return tri.Area();
		}else{ // either no intersection area, or circle in triangle
			if(inside & (1 << 3)){ // triangle contains circle center, intersection area is either circle area or triangle area
				return circle.Area();
			}else{
				return NumericType(0);
			}
		}
	}else if(2 == nx){
		// Either the 2 intersections are a single side or on two different sides
		if(nxp[0] < 2 && nxp[1] < 2 && nxp[2] < 2){ // on different sides
			// The area is determined by tracing
		}else{
			int i;
			for(i = 0; i < 3; ++i){
				if(nxp[i] > 1){ break; }
			}
			// Either the circle is mostly inside with a wedge poking out a side
			// or the circle is mostly outside with a wedge poking inside
			NumericType sector_area = CirclularSectorArea(circle.radius, (xp[2*i+1]-xp[2*i]).Length());
			if(inside & (1 << 3)){
				// Area of circle minus a wedge
				return circle.Area() - sector_area;
			}else{
				return sector_area;
			}
		}
	}else if(4 == nx){
		// The area is determined by tracing
	}else if(6 == nx){
		// The area is determined by tracing
	}else{
		// There is no way we can get here
		return NumericType(-1);
	}
	
	// At this point we expect to just trace out the intersection shape
	// The vertices of the intersection shape is either a triangle vertex
	// or a intersection point on a triangle edge.
	int vtype[6]; // 1 = triangle vertex, 0 = intersection point
	pt_t vp[6];
	int nv = 0; // number of actual vertices
	
	for(int i = 0; i < 3; ++i){
		if(inside & (1 << i)){
			vp[nv] = seg[i][0];
			vtype[nv++] = 1;
		}
		for(int j = 0; j < nxp[i]; ++j){
			vp[nv] = xp[2*i+j];
			vtype[nv++] = 0;
		}
	}
	
	if(nv < 3){ // this should not be possible
		return NumericType(-1);
	}
	
	// All neighboring points in v which are intersection points should have circular caps added
	NumericType area(0);
	for(int i = 2; i < nv; ++i){
		area += TTriangle2<NumericType>(vp[0], vp[i-1], vp[i]).Area();
		if((0 == vtype[i-1]) && (0 == vtype[i])){
			area += CirclularSectorArea(circle.radius, (vp[i]-vp[i-1]).Length());
		}
	}
	// Check the final segments (those next to vp[0]) to see if they need caps added
	if(0 == vtype[0]){
		if(0 == vtype[1]){
			area += CirclularSectorArea(circle.radius, (vp[1]-vp[0]).Length());
		}
		if(0 == vtype[nv-1]){
			area += CirclularSectorArea(circle.radius, (vp[0]-vp[nv-1]).Length());
		}
	}
	return area;
}

template <typename NumericType>
NumericType IntersectionArea(const TParallelogram2<NumericType> &quad, const TCircle2<NumericType> &circle){
	return
		IntersectionArea(TTriangle2<NumericType>(quad[0], quad[1], quad[2]), circle) +
		IntersectionArea(TTriangle2<NumericType>(quad[0], quad[1], quad[3]), circle);
}

template <typename NumericType>
NumericType IntersectionArea(const TTriangle2<NumericType> &T1, const TTriangle2<NumericType> &T2){
	typedef TSegment2<NumericType> segment_t;
	typedef TPt2<NumericType> pt_t;
	
	// This is much easier than the circle-triangle test; the insidedness of points is all that matters
	int inside = 0; // bit array of which T1 vertices are in T2
	for(i = 0; i < 3; ++i){
		if(T2.Contains(T1[i])){ inside |= (1 << i); }
	}
	
	// Get the 3 segments of the T1
	segment_t seg[3];
	for(size_t i = 0; i < 3; ++i){
		seg[i] = segment_t(T1[i], T2[i+1]);
	}
	
	// Get intersections of each segment of T1 with each T2
	pt_t xp[6]; // intersection points
	int nxp[3]; // number of intersections with each segment
	int nx = 0;
	for(size_t i = 0; i < 3; ++i){
		nxp[i] = Intersection(T2, seg[i], &(xp[2*i]));
		nx += nxp[i];
	}
	
	if(0 == nx){ // disjoint or one contained in other
		if(0x7 == inside){
			return T1.Area();
		}else if(T1.Contains(T2[0])){ // if T1 contains any vertex of T2, all of T2 is inside
			return T2.Area();
		}else{
			return NumericType(0);
		}
	}
	
	// Gather all interestion area vertices
	pt_t vp[6]; int nv = 0;
	for(int i = 0; i < 3; ++i){
		if(inside & (1 << i)){
			vp[nv++] = seg[i][0];
		}
		for(int j = 0; j < nxp[i]; ++j){
			vp[nv++] = xp[2*i+j];
		}
	}
	
	if(nv < 3){ return 0; }
	
	// All neighboring points in v which are intersection points should have circular caps added
	NumericType area(0);
	for(int i = 2; i < nv; ++i){
		area += TTriangle2<NumericType>(vp[0], vp[i-1], vp[i]).Area();
	}
	return area;
}

#endif // _TINTERSECTION2_HPP_
