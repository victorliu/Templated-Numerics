#ifndef _PERIODIC_TRIANGULATION2_HPP_
#define _PERIODIC_TRIANGULATION2_HPP_

#include <cstdlib>
#include <cmath>
#include <vector>
#include <AnalyticGeometry/TPt2.h>
#include <AnalyticGeometry/TVec2.h>

namespace Periodic{

template <typename NumericType>
class Triangulation{
public:
	class Vertex{
	public:
		size_t v;
		int offset[2];
	};
	class Edge{
	public:
		Vertex v[2];
		size_t equiv; // index of an equivalent edge
	};
	class Triangle{
	public:
		Edge e[3];
	};
private:
	Lattice2<NumericType> L;
public:
	Triangulation(const Lattice2<NumericType> &lattice, const std::vector<Pt2> &uv):L(lattice){
	}
};

}; // namespace Periodic

#endif // _PERIODIC_TRIANGULATION2_HPP_

