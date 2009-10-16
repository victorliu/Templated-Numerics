#ifndef _TINTERSECTION_HPP_
#define _TINTERSECTION_HPP_

// Returns number of intersections between circle and segment p0,p1, returns points in x0,x1
// x0 is closer than x1 to p0 if there are two intersections
template <typename NumericType>
int Intersect(const TCircle2<NumericType> &circle, const TSegment2<NumericType> &seg, TPt2<NumericType> *x[2]){
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
	if(d < 0){ return 0; }
	else if(d == 0){
		if(NULL != x0){
			x[0] = seg[0] - b_2a*p1p0;
		}
		return 1;
	}else{
		if(NULL != x0){
			x[0] = seg[0] + (-b_2a-sqrt(d))*p1p0;
		}
		if(NULL != x1){
			x[1] = seg[0] + (-b_2a+sqrt(d))*p1p0;
		}
		return 2;
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
	int last_xi = -1;
	for(size_t i = 0; i < 3; ++i){
		nxp[i] = Intersection(circle, seg[i], &(xp[2*i]));
		if(nxp[i] > 0){
			last_xi = i;
			nx += nxp[i];
		}
	}
	
	// The total number of intersections has to be even
	if(nx & 1){
		// Attempt to fix it; this can happen if an intersection is very near a vertex
		// Remove the intersection that is closest to a vertex
		int min_dist_xi = 2*tri.LongestSideLength();
		for(int i = 0; i < 3; ++i){
			for(int j = 0; j < nxp[i]; ++j){
			}
		}
	}
	
	int inside = 0; // bit array of which quad points are in circle, 5th bit is if circle center is in quad
	for(i = 0; i < 3; ++i){
		if(circle.Contains(quad[i])){ inside |= (1 << i); }
	}
	if(tri.Contains(circle.center)){ inside |= (1 << 3); }
	
	if(nx == 0){ // either no intersection area, or triangle entirely in circle, or circle in triangle
		if((inside & 0x7) == 0x7){ // all triangle points in circle
			return tri.Area();
		}else{ // either no intersection area, or circle in triangle
			if(inside & (1 << 3)){ // triangle contains circle center, intersection area is either circle area or triangle area
				return circle.Area();
			}else{
				return NumericType(0);
			}
		}
	}
	
	// At this point, the circle and triangle line shapes intersect
	// We need to "triangulate" the interior.
	// The interior is a polygon where each face might have a circular cap
	pt_t v[6]; int nv = 0;
	int vc[6] = {0,0,0,0,0,0}; // whether vc[i] = 1 if segment (v[i],v[i+1]) has a circular cap
	v[0] = xp[2*last_xi]; // this is a known intersection point from above
}

template <typename NumericType>
NumericType IntersectionArea(const TParallelogram2<NumericType> &quad, const TCircle2<NumericType> &circle){
	return
		IntersectionArea(TTriangle2<NumericType>(quad[0], quad[1], quad[2]), circle) +
		IntersectionArea(TTriangle2<NumericType>(quad[0], quad[1], quad[3]), circle);
}

#endif // _TINTERSECTION_HPP_
