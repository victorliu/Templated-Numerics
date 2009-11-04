#define _USE_MATH_DEFINES
#include <iostream>
//#define USE_INTEGRATION_FIXED_DIMENSION
#include "Integrator_meta.h"

//#define which 0
//#define dim 3
const double radius = 0.50124145262344534123412; /* random */

struct MyType{
	typedef double value_type;
	double val;
	double payload;
	MyType(double v = 0):val(v),payload(0){}
	MyType(double v, double p):val(v),payload(p){}
	MyType& operator+=(const MyType &other){
		val += other.val;
		payload += other.payload;
		return *this;
	}
	MyType& operator-=(const MyType &other){
		val -= other.val;
		//payload -= other.payload;
		return *this;
	}
	MyType& operator*=(const value_type &scale){
		val *= scale;
		//payload *= scale;
		return *this;
	}
	operator double() const{
		return val;
	}
};
MyType operator*(const MyType::value_type &a, const MyType &b){
	return (MyType(b) *= a);
}
MyType operator+(const MyType &a, const MyType &b){
	return (MyType(a) += b);
}

struct MyFunc :
#ifdef USE_INTEGRATION_FIXED_DIMENSION
	public Integration::FunctionN<dim,double,MyType>
#else
	public Integration::Function<double,MyType>
#endif
{
#ifdef USE_INTEGRATION_FIXED_DIMENSION
	MyType operator()(const double *x) const
#else
	MyType operator()(const std::vector<double> &x) const
#endif
	{
		double val;
		int i;
		double scale = 1.0;

		switch(which){
		case 0:
			val = 0;
			for (i = 0; i < dim; ++i)
				val += x[i] * x[i];
			val = (val < radius * radius) ? 1 : 0;
			break;
		case 1:
			val = 1;
			for (i = 0; i < 2; ++i) val *= cos(x[i]);
			break;
		case 2:
			val = 0;
			for (i = 0; i < dim; ++i) {
				double z = (1 - x[i]) / x[i];
				val += z * z;
				scale *= M_2_SQRTPI / (x[i] * x[i]);
			}
			val = exp(-val) * scale;
			break;
		default:
			val = 1.0-(x[0]*x[0]+x[1]*x[1]);
			if(val < 0){ val = 0; }
			std::cout << x[0] << "\t" << x[1] << "\t" << val << std::endl;
			break;
		}
		return MyType(val,1.0);
	}
	size_t Dimension() const{ return dim; }
	MyFunc(size_t n, size_t which_integrand):dim(n),which(which_integrand){}
private:
	size_t dim, which;
};

#ifdef USE_INTEGRATION_FIXED_DIMENSION
typedef Integration::Cubature<dim,MyFunc> Integrator2;
#else
typedef Integration::Cubature<MyFunc> Integrator2;
#endif

int main(int argc, char **argv){
     size_t ndim = argc > 1 ? atoi(argv[1]) : 2;
     double tol = argc > 2 ? atof(argv[2]) : 1e-2;
     size_t which_integrand = argc > 3 ? atoi(argv[3]) : 1;
     size_t maxEval = argc > 4 ? atoi(argv[4]) : 0;


	Integrator2::Parameters params;
	params.max_evals = maxEval;
	params.max_absolute_error = 0;
	params.max_relative_error = tol;
	MyFunc::Domain dom(ndim);
	for(size_t i = 0; i < ndim; ++i){
		dom.a[i] = 0;
		if(which_integrand == 2){
			dom.b[i] = 1;
		}else{
			dom.b[i] = 1 + 0.4*sin((double)i);
		}
	}
	
	Integrator2 integrator(MyFunc(ndim,which_integrand), dom, params);
	double error;
	MyType result;
	integrator.Integrate(result, &error);
	
	std::cout << "result = " << result << ", error = " << error << std::endl;
	
	return 0;
}

