#ifndef _TPREDICATES2_HPP_
#define _TPREDICATES2_HPP_

// Preprocessor flags:
//   USE_SHEWCHUK_PREDICATES (otherwise use the fast approximate method)
//     Uses J. R. Shewchuk's adaptive precision primitives (public domain)
//   USE_SHEWCHUK_PREDICATES_TYPE (defaults to double)
//   USE_ABDPY (cannot use with USE_SHEWCHUK_PREDICATES)
//     Uses exact determinant sign computation by F. Avnaim, et al.
//     Not integrated yet

#include "TSegment2.h"
#include "TTriangle2.h"
#include "TVec2.h"

// Standing at segment[0], looking towards segment[1], is p on the left?
// (positive for left, negative for right, zero for exactly on)
template <typename NumericType>
NumericType Orient2(const TSegment2<NumericType> &segment, const TPt2<NumericType> &p){
	return TVec2<NumericType>::Cross(segment[1]-segment[0], p-segment[0]);
}
template <typename NumericType>
NumericType Orient2(const TPt2<NumericType> &a, const TPt2<NumericType> &b, const TPt2<NumericType> &p){
	return TVec2<NumericType>::Cross(b-a, p-a);
}
template <typename NumericType>
NumericType Orient2(const TLine2<NumericType> &line, const TPt2<NumericType> &p){
	return TVec2<NumericType>::Cross(line.v, p-line.p);
}

// Is p in the circumcircle of (a,b,c)? (positive for yes, negative for outside, zero is on)
template <typename NumericType>
NumericType InCircle2(const TPt2<NumericType> &a,
                      const TPt2<NumericType> &b,
                      const TPt2<NumericType> &c,
                      const TPt2<NumericType> &p)
{
	typedef TVec2<NumericType> Vec2;
	Vec2 ap(a-p);
	Vec2 bp(b-p);
	Vec2 cp(c-p);
	
	NumericType abdet(Vec2::Cross(ap, bp));
	NumericType bcdet(Vec2::Cross(bp, cp));
	NumericType cadet(Vec2::Cross(cp, ap));
	
	return bcdet*ap.LengthSq() + cadet*bp.LengthSq() + abdet*cp.LengthSq();
}
template <typename NumericType>
NumericType InCircle2(const TTriangle2<NumericType> &t,
                      const TPt2<NumericType> &p)
{
	return InCircle2(t[0], t[1], t[2], p);
}

#ifdef USE_SHEWCHUK_PREDICATES
# ifndef USE_SHEWCHUK_PREDICATES_TYPE
#  define USE_SHEWCHUK_PREDICATES_TYPE double
# endif

extern "C" USE_SHEWCHUK_PREDICATES_TYPE orient2d(
	USE_SHEWCHUK_PREDICATES_TYPE* pa,
	USE_SHEWCHUK_PREDICATES_TYPE* pb,
	USE_SHEWCHUK_PREDICATES_TYPE* pc
	);
extern "C" USE_SHEWCHUK_PREDICATES_TYPE incircle(
	USE_SHEWCHUK_PREDICATES_TYPE* pa,
	USE_SHEWCHUK_PREDICATES_TYPE* pb,
	USE_SHEWCHUK_PREDICATES_TYPE* pc,
	USE_SHEWCHUK_PREDICATES_TYPE* pd
	);

template <>
USE_SHEWCHUK_PREDICATES_TYPE Orient2(
	const TSegment2<USE_SHEWCHUK_PREDICATES_TYPE> &segment,
	const TPt2<USE_SHEWCHUK_PREDICATES_TYPE> &p)
{
	return orient2d(
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&segment[0]),
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&segment[1]),
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&p)
		);
}
template <>
USE_SHEWCHUK_PREDICATES_TYPE Orient2(
	const TPt2<USE_SHEWCHUK_PREDICATES_TYPE> &a,
	const TPt2<USE_SHEWCHUK_PREDICATES_TYPE> &b,
	const TPt2<USE_SHEWCHUK_PREDICATES_TYPE> &p)
{
	return orient2d(
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&a),
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&b),
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&p)
		);
}
template <>
USE_SHEWCHUK_PREDICATES_TYPE Orient2(
	const TLine2<USE_SHEWCHUK_PREDICATES_TYPE> &line,
	const TPt2<double> &p)
{
	return Orient2(TSegment2<USE_SHEWCHUK_PREDICATES_TYPE>(line.p, line.p+line.v), p);
}
template <>
USE_SHEWCHUK_PREDICATES_TYPE InCircle2(
	const TPt2<USE_SHEWCHUK_PREDICATES_TYPE> &a,
	const TPt2<USE_SHEWCHUK_PREDICATES_TYPE> &b,
	const TPt2<USE_SHEWCHUK_PREDICATES_TYPE> &c,
	const TPt2<USE_SHEWCHUK_PREDICATES_TYPE> &p)
{
	return incircle(
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&a),
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&b),
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&c),
		reinterpret_cast<USE_SHEWCHUK_PREDICATES_TYPE*>(&p)
		);
}
#endif // USE_SHEWCHUK_PREDICATES

#endif // _TPREDICATES2_HPP_
