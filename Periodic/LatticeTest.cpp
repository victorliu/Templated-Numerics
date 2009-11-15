#include <NumericTypes/rational.h>
#include <NumericTypes/rational_radical1.hpp>
#include "Lattice2.hpp"

typedef rational_radical1<3> real_type;
typedef Lattice2<real_type> Lattice;
typedef Lattice::Vec2 Vec2;

void print_out(const Lattice &L){
	Lattice::LatticeSymmetry type;
	std::vector<Lattice::Mat2> symmetries;
	
	type = L.GetSymmetries(symmetries);
	
	Lattice::Vec2 u,v;
	L.GetBasis(u,v);
	std::cout << "Basis: " << u << " " << v << std::endl;
	
	std::cout << "Lattice type: ";
	switch(type){
	case Lattice::OBLIQUE:
		std::cout << "Oblique" << std::endl;
		break;
	case Lattice::RHOMBIC:
		std::cout << "Rhombic" << std::endl;
		break;
	case Lattice::RECTANGULAR:
		std::cout << "Rectangular" << std::endl;
		break;
	case Lattice::SQUARE:
		std::cout << "Square" << std::endl;
		break;
	case Lattice::HEXAGONAL:
		std::cout << "Hexagonal" << std::endl;
		break;
	}
	std::cout << std::endl;
	
	std::cout << "Lattice symmetry transformations:" << std::endl;
	for(size_t i = 0; i < symmetries.size(); ++i){
		std::cout << symmetries[i] << std::endl;
	}
	std::cout << std::endl;
	
	std::cout << "Ring iteration test:" << std::endl;
	int a = 0, b = 0;
	std::vector<Lattice::Coord2> points;
	for(int i = 0; i < 15; ++i){
		points.clear();
		L.GetNextPointRing(a, b, points);
		for(int j = 0; j < points.size(); ++j){
			std::cout << points[j] << " ";
		}
		std::cout << std::endl;
	}
	
	std::cout << "MP Special points test:" << std::endl;
	std::vector<Lattice::SpecialPoint> spoints;
	for(int i = 1; i <= 19; i+=2){
		L.GetMonkhorstPackPoints(i, spoints);
		for(int j = 0; j < spoints.size(); ++j){
			//std::cout << "{" << spoints[j].x << "," << spoints[j].weight << "} ";
			std::cout << spoints[j].x << ", ";
		}
		std::cout << std::endl;
	}
	
	std::cout << "---" << std::endl;
}

int main(){
/*
Oblique:        (3,0) (1,2)
Square:         (1,0) (0,1)               u.v = 0, |u| = |v|
Rectangular:    (1,0) (0,2)               u.v = 0
Hexagonal:      (1,0) (1/2,sqrt(3)/2)    |u.v| = 1/2, |u| = |v|
Centered rect:  (1,0) (1/2, 1)           |u.v| = 1/2 or equivalently |u| = |v|
*/
	
	rational half(1,2);
	
	Lattice l01(Vec2(3,0), Vec2(1,2));
	print_out(l01);
	
	Lattice l02(Vec2(0,1), Vec2(1,0));
	print_out(l02);
	
	Lattice l03(Vec2(0,2), Vec2(1,0));
	print_out(l03);
	
	Lattice l04(Vec2(1,0), Vec2(half,real_type(0,half)));
	print_out(l04);
	
	Lattice l05(Vec2(2,0), Vec2(1,2));
	print_out(l05);
	
	Lattice l06(Vec2(2,1), Vec2(1,2));
	print_out(l06);
	
	Lattice l07(Vec2(1,98), Vec2(1,99));
	print_out(l07);
	
	Lattice l08(Vec2(2,real_type(0,1)), Vec2(1,real_type(0,1)));
	print_out(l08);

	return 0;
}
