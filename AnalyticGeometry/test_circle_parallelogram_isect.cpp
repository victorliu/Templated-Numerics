#include <iostream>
#include <AnalyticGeometry/TCircle2.h>
#include <AnalyticGeometry/TParallelogram2.h>
#include <AnalyticGeometry/TIntersection2.hpp>

typedef double real_t;
typedef TPt2<real_t> pt_t;
typedef TVec2<real_t> vec_t;
typedef TCircle2<real_t> circle_t;
typedef TParallelogram2<real_t> parallelogram_t;

int main(){
	circle_t C(pt_t::Origin, 5);
	for(int i = -10; i <= 10; ++i){
	//{int i = -6;
		for(int j = -10; j <= 10; ++j){
		//{int j = -1;
			parallelogram_t cell(pt_t(real_t(i),real_t(j)), vec_t(1,0), vec_t(0,1));
			real_t area = IntersectionArea(cell, C);
			std::cout << i << "\t" << j << "\t" << area << std::endl;
		}std::cout << std::endl;
	}
	return 0;
}
