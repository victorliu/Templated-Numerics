#ifndef _INTERPOLATE_HPP_
#define _INTERPOLATE_HPP_

#include <AnalyticGeometry/TVec2.hpp>
#include <Periodic/Lattice2.hpp>

namespace Interpolation{

enum Method{
	LINEAR,
	QUADRATIC,
	CUBIC,
	SINC
};

template <typename LatticeType, typename NumericType>
NumericType Interpolate_Linear(const Lattice2<LatticeType> &L, const NumericType *grid, size_t n1, size_t n2,
                 TVec2<LatticeType> &uvcoords)
{
	// Bring uvcoords within the parallelogram defined by the grid
	uvcoords[0] = uvcoords[0] - floor(uvcoords[0]);
	uvcoords[1] = uvcoords[1] - floor(uvcoords[1]);
	
	
}

// Periodic interpolation over a grid on a lattice
// L is the lattice, so L.u is divided into n1 pieces, and L.v is divided into n2 pieces
// grid is the values over the 2D point array over the fundamental parallelogram of L
// uvcoords is the point in the coordinates of the L frame in which we want to interpolate
template <typename LatticeType, typename NumericType>
NumericType Interpolate(const Lattice2<LatticeType> &L, const NumericType *grid, size_t n1, size_t n2,
                 TVec2<LatticeType> &uvcoords, Method method = LINEAR)
{
	switch(method){
	case LINEAR:
		return Interpolate_Linear(L, grid, n1, n2, uvcoords);
	default:
		return Interpolate_Linear(L, grid, n1, n2, uvcoords);
	}
}

}; // namespace Interpolation

#endif // _INTERPOLATE_HPP_
