#ifndef _TGRID2_H_
#define _TGRID2_H_

#include "TArray2.h"
#include <AnalyticGeometry/TVec2.h>
#include <AnalyticGeometry/TPt2.h>

template <typename CoordType, typename ValueType, class TAllocator = std::allocator<ValueType> >
class TGrid2{
public:
	typedef ValueType value_type;
	typedef TArray2<value_type,TAllocator> Array2;
	typedef TVec2<CoordType> Vec2;
	typedef TPt2<CoordType> Pt2;
	
	TGrid2(const Pt2 &Base, const Vec2 &Vec0, const Vec2 &Vec1, size_t m, size_t n, const value_type& init_val = value_type()):A(m,n,init_val),u(Vec0),v(Vec1),base(Base){
	}
	TGrid2(const TGrid2 &g):A(g.A),u(g.u),v(g.v),base(g.base){
	}
	TGrid2& operator=(const TGrid2 &g){
		if(this != &g){
			A = g.A;
			u = g.u;
			v = g.v;
			base = g.base;
		}
		return *this;
	}
	
	void Resize(size_t N0, size_t N1){
		A.Resize(N0, N1);
	}
	
	// GetPoint does not check to see if the coordinates (i,j) are in range
	// This allows generation of any grid point
	inline Pt2 GetPoint(int i, int j) const{ return base + (CoordType(i)/CoordType(Dim0()-1))*u + (CoordType(j)/CoordType(Dim1()-1))*v; }
	inline const value_type& operator()(size_t i, size_t j) const{ return A(i,j); }
	inline       value_type& operator()(size_t i, size_t j)      { return A(i,j); }
	inline const Vec2& Side0() const{ return u; }
	inline const Vec2& Side1() const{ return v; }
	inline Vec2 Inc0() const{ return u/CoordType(Dim0()); }
	inline Vec2 Inc1() const{ return v/CoordType(Dim1()); }
	inline size_t Dim0() const{ return A.Dim0(); }
	inline size_t Dim1() const{ return A.Dim1(); }
	inline value_type* Raw(){ return A.Raw(); }
	inline const Array2& Array() const{ return A; }
private:
	TArray2<value_type,TAllocator> A;
	Vec2 u, v;
	Pt2 base;
};

#endif // _TGRID2_H_
