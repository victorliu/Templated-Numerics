#include "rational.h"
#include "rational_radical1.hpp"
#include "TVec2.h"
#include "TMat2.h"
#include <set>

typedef rational_radical1<7> rr3;
typedef TVec2<rr3> tvec;
typedef TVec2<float> fvec;

float xform(float x){
	return 10*x;
}
fvec xformp(const tvec &r){
	return fvec(
		300+xform(r.v[0].float_value()),
		400+xform(r.v[1].float_value())
		);
}

void draw_point(const tvec &r){
	const float rad = 1;
	std::cout << "newpath" << std::endl;
	fvec p0, p1;
	p0 = xformp(r); p1 = p0;
	p0.v[0] -= rad; p1.v[0] += rad;
	std::cout << p0.v[0] << ' ' << p0.v[1] << " moveto" << std::endl;
	std::cout << p1.v[0] << ' ' << p1.v[1] << " lineto" << std::endl;
	std::cout << "stroke" << std::endl;
	
	std::cout << "newpath" << std::endl;
	p0 = xformp(r); p1 = p0;
	p0.v[1] -= rad; p1.v[1] += rad;
	std::cout << p0.v[0] << ' ' << p0.v[1] << " moveto" << std::endl;
	std::cout << p1.v[0] << ' ' << p1.v[1] << " lineto" << std::endl;
	std::cout << "stroke" << std::endl;
}

void draw_circle(float r){
	fvec center = xformp(tvec(0,0));
	std::cout << center.v[0] << ' ' << center.v[1] << ' ' << xform(r) << " 0 360 arc closepath" << std::endl;
	std::cout << "stroke" << std::endl;
}

int main(){
	tvec u(1,0), v(0,1);
	std::set<rr3> r2set;
	
	for(int j = -3; j <= 3; ++j){
		for(int i = -6; i <= 6; ++i){
			tvec x = rr3(i)*u + rr3(j)*v;
			rr3 l2 = x.LengthSq();
			draw_point(x);
			r2set.insert(l2);
			//std::cout << l2.float_value() << "\t";
		}
		//std::cout << std::endl;
	}
	//*
	for(std::set<rr3>::const_iterator i = r2set.begin(); i != r2set.end(); ++i){
		draw_circle(sqrt(i->float_value()));
	}
	//*/
	
	return 0;
}
