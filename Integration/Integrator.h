#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

// Licensing details:
//   The cubature code is under GNU General Public License

#include <vector>
#include <NumericTypes/Traits.hpp>
#include <queue>

// Preprocessor macros
//   USE_INTEGRATION_STATUS_CALLBACK
//   USE_INTEGRATION_FIXED_DIMENSION

// If DomainType in general must be a builtin type (double or float)
// since it must have an available std::numeric_traits<>.
// ReturnType can be any type for which the following conditions hold:
//   Must have a default constructor which initializes the DomainType cast value to zero
//   Defines operator DomainType() const - the type conversion operator to extract the value to be resolved
//   Defines operator+=, operator-=, and operator*=
//   There exists operator+(const T&, const T&) and operator*(const DomainType&, const T&)

#ifdef USE_INTEGRATION_STATUS_CALLBACK
# include "Status.h"
#endif

namespace Integration{

template <unsigned int N, class DomainType, class ReturnType>
class FunctionN{
public:
	typedef DomainType value_type;
	typedef ReturnType return_type;
	struct Domain{
		value_type a[N], b[N];
	};
	virtual return_type operator()(const value_type *x) const = 0;
	size_t Dimension() const{ return N; }
	virtual size_t MaxParallel() const{ return 1; }
};
template <class DomainType, class ReturnType>
class Function{
public:
	typedef DomainType value_type;
	typedef ReturnType return_type;
	struct Domain{
		std::vector<value_type> a, b;
		Domain(size_t dim):a(dim),b(dim){}
	};
	virtual return_type operator()(const std::vector<value_type> &x) const = 0;
	virtual size_t Dimension() const = 0;
	virtual size_t MaxParallel() const{ return 1; }
};

template <unsigned int N, typename FunctionType>
class IntegratorN{
public:
	typedef FunctionType function_type;
	typedef typename FunctionType::Domain domain_type;
	typedef typename FunctionType::value_type value_type;
	typedef typename FunctionType::return_type return_type;
	// ReturnType must have an operator value_type() conversion function;
	
	
	IntegratorN(function_type func, const domain_type &dom
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		, const StatusCallback &status_callback = DefaultStatusCallback()
#endif
		):
		f(func),d(dom)
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		, status_func(status_callback)
#endif
	{
	}
	virtual ~IntegratorN(){}
	virtual int Integrate(return_type &val, value_type *error = NULL) = 0;
	virtual int Refine(return_type &val, value_type *error = NULL) = 0;
protected:
	function_type f;
	domain_type d;
#ifdef USE_INTEGRATION_STATUS_CALLBACK
	StatusCallback status_func;
#endif
};

template <typename FunctionType>
class Integrator{
public:
	typedef FunctionType function_type;
	typedef typename FunctionType::Domain domain_type;
	typedef typename FunctionType::value_type value_type;
	typedef typename FunctionType::return_type return_type;
	
	Integrator(function_type func, const domain_type &dom
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		, const StatusCallback &status_callback = DefaultStatusCallback()
#endif
		):
		f(func), d(dom)
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		, status_func(status_callback)
#endif
	{
	}
	virtual ~Integrator(){}
	virtual int Integrate(return_type & val, value_type *error = NULL) = 0;
	virtual int Refine(return_type &val, value_type *error = NULL) = 0;
protected:
	function_type f;
	domain_type d;
#ifdef USE_INTEGRATION_STATUS_CALLBACK
	StatusCallback status_func;
#endif
};

template <
#ifdef USE_INTEGRATION_FIXED_DIMENSION
	unsigned int N,
#endif
	typename FunctionType>
class Cubature :
#ifdef USE_INTEGRATION_FIXED_DIMENSION
	public IntegratorN<N,FunctionType>
#else
	public Integrator<FunctionType>
#endif
{
public:
	typedef FunctionType function_type;
	typedef typename FunctionType::Domain domain_type;
	typedef typename FunctionType::value_type value_type;
	typedef typename FunctionType::return_type return_type;
	
	struct Parameters{
		size_t max_evals;
		value_type max_absolute_error;
		value_type max_relative_error;
	};
	
	Cubature(
		function_type func, const domain_type &dom,
		const Parameters &params
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		, const StatusCallback &status_callback = DefaultStatusCallback()
#endif
		):
#ifdef USE_INTEGRATION_FIXED_DIMENSION
		IntegratorN<N,function_type>
#else
		Integrator<function_type>
#endif
			(func, dom
#ifdef USE_INTEGRATION_STATUS_CALLBACK
			, status_callback
#endif
		),
		_params(params),h(dom.a, dom.b)
	{
		const size_t n = func.Dimension();
		h.volume = 1;
		for(size_t i = 0; i < n; ++i){
			h.center[i] = (dom.a[i]+dom.b[i])/value_type(2);
			h.half_width[i] = (dom.b[i]-dom.a[i])/value_type(2);
			h.volume *= 2*h.half_width[i];
		}
	}
	int Integrate(return_type &val, value_type *error = NULL){
		const size_t dim = this->f.Dimension();
		Rule *r;
		switch(dim){
		case 1:
			r = new Rule15Gauss();
			break;
		default:
			r = new Rule75GenzMalik(dim);
			break;
		}
		int status = ruleadapt_integrate(r, h, val, error);
		delete r;
		return status;
	}
	int Refine(return_type &val, value_type *error = NULL){
		return 0;
	}
private:
	struct ValErr{
		value_type err;
		return_type val;
		ValErr():err(0){}
	};
	struct Hypercube{
#ifdef USE_INTEGRATION_FIXED_DIMENSION
		value_type center[N];
		value_type half_width[N];
#else
		std::vector<value_type> center, half_width;
#endif
		value_type volume;
#ifdef USE_INTEGRATION_FIXED_DIMENSION
		Hypercube(const value_type *centers, const value_type *hwidths)
#else
		Hypercube(const std::vector<value_type> &centers, const std::vector<value_type> &hwidths):center(centers),half_width(hwidths)
#endif
		{
#ifdef USE_INTEGRATION_FIXED_DIMENSION
			const size_t n = N;
#else
			const size_t n = centers.size();
#endif
			volume = value_type(1);
			for(size_t i = 0; i < n; ++i){
				volume *= value_type(2)*half_width[i];
#ifdef USE_INTEGRATION_FIXED_DIMENSION
				center[i] = centers[i];
				half_width[i] = hwidths[i];
#endif
			}
		}
		Hypercube& operator=(const Hypercube &h){
#ifdef USE_INTEGRATION_FIXED_DIMENSION
			const size_t n = N;
#else
			const size_t n = h.center.size();
#endif
			for(size_t i = 0; i < n; ++i){
				center[i] = h.center[i];
				half_width[i] = h.half_width[i];
			}
			volume = h.volume;
			return *this;
		}
	};
	struct Rule;
	struct Region{
		Hypercube h;
		ValErr ve;
		size_t split_dim;
		Region(const Hypercube &cube):h(cube),split_dim(0){}
		void Split(Region &new_region){
			size_t d = split_dim;
			new_region.ve = ve;
			new_region.split_dim = split_dim;
			h.half_width[d] /= return_type(2);
			h.volume /= return_type(2);
			new_region.h = h;
			h.center[d] -= h.half_width[d];
			new_region.h.center[d] += h.half_width[d];
		}
		void Eval(const function_type &f, const Rule &r){
			split_dim = r.EvalError(f, h, ve);
		}
		bool operator<(const Region &b) const{
			return (ve.err < b.ve.err);
		}
	};
	struct Rule{
	public:
		size_t dim;
		size_t num_points;
		Rule(size_t Dim, size_t npoints = 0):dim(Dim),num_points(npoints){}
		virtual size_t EvalError(const function_type &f, const Hypercube &h, ValErr &ve) const = 0;
	};
	// ls0 returns the least-significant 0 bit of n (e.g. it returns
	// 0 if the LSB is 0, it returns 1 if the 2 LSBs are 01, etc.).
	static size_t ls0(size_t n){
		const unsigned bits[] = {
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
			0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8,
		};
		size_t bit;
		bit = 0;
		while((n & 0xff) == 0xff){
			n >>= 8;
			bit += 8;
		}
		return bit + bits[n & 0xff];
	}
	static return_type evalR_Rfs(const function_type &f,
#ifdef USE_INTEGRATION_FIXED_DIMENSION
	value_type *cur_point, const value_type *c, const value_type *r
#else
	std::vector<value_type> &cur_point, const std::vector<value_type> &c, const std::vector<value_type> &r
#endif
	){
		return_type sum(0);
		size_t signs = 0;	/* 0/1 bit = +/- for corresponding element of r[] */
		const size_t n = f.Dimension();

		/* We start with the point where r is ADDed in every coordinate (This implies signs=0) */
		for(size_t i = 0; i < n; ++i){
			cur_point[i] = c[i] + r[i];
		}

		/* Loop through the points in gray-code ordering */
		for(size_t i = 0;; ++i){
			size_t mask;

			sum += f(cur_point);

			size_t d = ls0(i);	/* which coordinate to flip */
			if(d >= n){ break; }

			/* flip the d-th bit and add/subtract r[d] */
			mask = 1U << d;
			signs ^= mask;
			cur_point[d] = (signs & mask) ? c[d] - r[d] : c[d] + r[d];
		}
		return sum;
	}
	
	static return_type evalRR0_0fs(const function_type &f,
#ifdef USE_INTEGRATION_FIXED_DIMENSION
	value_type *cur_point, const value_type *c, const value_type *r
#else
	std::vector<value_type> &cur_point, const std::vector<value_type> &c, const std::vector<value_type> &r
#endif
	){
		return_type sum = 0;
		const size_t n = f.Dimension();

		for(size_t i = 0; i < n-1; ++i){
			cur_point[i] = c[i] - r[i];
			for(size_t j = i + 1; j < n; ++j){
				cur_point[j] = c[j] - r[j];
				sum += f(cur_point);
				cur_point[i] = c[i] + r[i];
				sum += f(cur_point);
				cur_point[j] = c[j] + r[j];
				sum += f(cur_point);
				cur_point[i] = c[i] - r[i];
				sum += f(cur_point);

				cur_point[j] = c[j];	// Done with j -> Restore p[j]
			}
			cur_point[i] = c[i];		// Done with i -> Restore p[i]
		}
		return sum;
	}

	static size_t evalR0_0fs4d(const function_type &f,
#ifdef USE_INTEGRATION_FIXED_DIMENSION
	value_type *cur_point, const value_type *c, const value_type *r1, const value_type *r2
#else
	std::vector<value_type> &cur_point, const std::vector<value_type> &c, const std::vector<value_type> &r1, const std::vector<value_type> &r2
#endif
	, return_type &sum0, return_type &sum1, return_type &sum2
	){
		size_t dimDiffMax = 0;
		const size_t n = f.Dimension();

		value_type ratio = r1[0] / r2[0];

		ratio *= ratio;
		sum0 = f(cur_point);

		value_type maxdiff = 0;
		for(size_t i = 0; i < n; i++){
			return_type f1a, f1b, f2a, f2b;

			cur_point[i] = c[i] - r1[i];
			sum1 += (f1a = f(cur_point));
			cur_point[i] = c[i] + r1[i];
			sum1 += (f1b = f(cur_point));
			cur_point[i] = c[i] - r2[i];
			sum2 += (f2a = f(cur_point));
			cur_point[i] = c[i] + r2[i];
			sum2 += (f2b = f(cur_point));
			cur_point[i] = c[i];

			value_type diff = std::abs(value_type(f1a) + value_type(f1b) - 2 * value_type(sum0) - ratio * (value_type(f2a) + value_type(f2b) - 2 * value_type(sum0)));

			if(diff > maxdiff){
				maxdiff = diff;
				dimDiffMax = i;
			}
		}

		return dimDiffMax;
	}

	struct Rule75GenzMalik : public Rule{
		// dimension-dependent constants
		value_type weight1, weight3, weight5;
		value_type weightE1, weightE3;
		 
		Rule75GenzMalik(size_t Dim):Rule(Dim, 1 + 2 * 2*Dim + 2*Dim*(Dim-1) + (1<<Dim)){
			// Needs dim >= 2
			weight1 = (value_type(12824 - 9120 * int(Dim) + 400 * int(Dim*Dim)) / value_type(19683));
			weight3 = value_type(1820 - 400 * int(Dim)) / value_type(19683);
			weight5 = value_type(6859) / value_type(19683) / value_type(1U << Dim);
			weightE1 = (value_type(729 - 950 * int(Dim) + 50 * int(Dim*Dim)) / value_type(729));
			weightE3 = value_type(265 - 100 * int(Dim)) / value_type(1458);
		}
		size_t EvalError(const function_type &f, const Hypercube &h, ValErr &ve) const{
			// lambda2 = sqrt(9/70), lambda4 = sqrt(9/10), lambda5 = sqrt(9/19)
			const value_type lambda2 = 0.3585685828003180919906451539079374954541;
			const value_type lambda4 = 0.9486832980505137995996680633298155601160;
			const value_type lambda5 = 0.6882472016116852977216287342936235251269;
			const value_type weight2 = 980. / 6561.;
			const value_type weight4 = 200. / 19683.;
			const value_type weightE2 = 245. / 486.;
			const value_type weightE4 = 25. / 729.;

			return_type sum1(0), sum2(0), sum3(0), result;
#ifdef USE_INTEGRATION_FIXED_DIMENSION
			value_type p[N];
			value_type widthLambda[N];
			value_type widthLambda2[N];
#else
			std::vector<value_type> p(f.Dimension());
			std::vector<value_type> widthLambda(f.Dimension());
			std::vector<value_type> widthLambda2(f.Dimension());
#endif

			for(size_t i = 0; i < f.Dimension(); ++i){
				p[i] = h.center[i];
			}

			for(size_t i = 0; i < f.Dimension(); ++i){
				widthLambda2[i] = h.half_width[i] * lambda2;
			}
			for(size_t i = 0; i < f.Dimension(); ++i){
				widthLambda[i] = h.half_width[i] * lambda4;
			}

			/* Evaluate function in the center, in f(lambda2,0,...,0) and
			f(lambda3=lambda4, 0,...,0).  Estimate dimension with largest error */
			size_t dimDiffMax = evalR0_0fs4d(f, p, h.center, widthLambda2, widthLambda, sum1, sum2, sum3);

			/* Calculate sum4 for f(lambda4, lambda4, 0, ...,0) */
			return_type sum4(evalRR0_0fs(f, p, h.center, widthLambda));

			/* Calculate sum5 for f(lambda5, lambda5, ..., lambda5) */
			for(size_t i = 0; i < f.Dimension(); ++i){
				widthLambda[i] = h.half_width[i] * lambda5;
			}
			return_type sum5(evalR_Rfs(f, p, h.center, widthLambda));

			/* Calculate fifth and seventh order results */

			result = h.volume * (weight1 * sum1 + weight2 * sum2 + weight3 * sum3 + weight4 * sum4 + weight5 * sum5);
			value_type res5th = h.volume * (weightE1 * value_type(sum1) + weightE2 * value_type(sum2) + weightE3 * value_type(sum3) + weightE4 * value_type(sum4));

			ve.val = result;
			ve.err = std::abs(res5th - value_type(result));

			return dimDiffMax;
		}
	};
	struct Rule15Gauss : public Rule{
		Rule15Gauss():Rule(1, 15){
			// Needs Dim == 1
		}
		size_t EvalError(const function_type &f, const Hypercube &h, ValErr &ve) const{
			// Gauss quadrature weights and kronrod quadrature abscissae and
			// weights as evaluated with 80 decimal digit arithmetic by
			// L. W. Fullerton, Bell Labs, Nov. 1981.
			const size_t n = 8;
			const value_type xgk[8] = {  // abscissae of the 15-point kronrod rule
				0.991455371120812639206854697526329,
				0.949107912342758524526189684047851,
				0.864864423359769072789712788640926,
				0.741531185599394439863864773280788,
				0.586087235467691130294144838258730,
				0.405845151377397166906606412076961,
				0.207784955007898467600689403773245,
				0.000000000000000000000000000000000
			// xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
			// xgk[0], xgk[2], ... to optimally extend the 7-point gauss rule
			};
			static const value_type wg[4] = {  // weights of the 7-point gauss rule
				0.129484966168869693270611432679082,
				0.279705391489276667901467771423780,
				0.381830050505118944950369775488975,
				0.417959183673469387755102040816327
			};
			static const value_type wgk[8] = { // weights of the 15-point kronrod rule
				0.022935322010529224963732008058970,
				0.063092092629978553290700663189204,
				0.104790010322250183839876322541518,
				0.140653259715525918745189590510238,
				0.169004726639267902826583426598550,
				0.190350578064785409913256402421014,
				0.204432940075298892414161999234649,
				0.209482141084727828012999174891714
			};
#ifdef USE_INTEGRATION_FIXED_DIMENSION
			value_type x[1];
#else
			std::vector<value_type> x(1);
#endif
			const value_type center = h.center[0];
			const value_type width = h.half_width[0];
			return_type fv1[n - 1], fv2[n - 1];
			x[0] = center;
			const return_type f_center = f(x);
			return_type result_gauss = wg[n/2 - 1] * f_center;
			return_type result_kronrod = wgk[n - 1] * f_center;
			value_type result_abs = fabs(value_type(result_kronrod));
			value_type result_asc, mean, err;
			unsigned j;

			for (j = 0; j < (n - 1) / 2; ++j) {
				int j2 = 2*j + 1;
				return_type f1, f2, fsum;
				value_type w = width * xgk[j2];
				x[0] = center - w; fv1[j2] = f1 = f(x);
				x[0] = center + w; fv2[j2] = f2 = f(x);
				fsum = f1 + f2;
				result_gauss += wg[j] * fsum;
				result_kronrod += wgk[j2] * fsum;
				result_abs += wgk[j2] * (abs(f1) + abs(f2));
			}

			for (j = 0; j < n/2; ++j) {
				int j2 = 2*j;
				return_type f1, f2;
				value_type w = width * xgk[j2];
				x[0] = center - w; fv1[j2] = f1 = f(x);
				x[0] = center + w; fv2[j2] = f2 = f(x);
				result_kronrod += wgk[j2] * (f1 + f2);
				result_abs += wgk[j2] * (abs(f1) + abs(f2));
			}

			ve.val = result_kronrod * width;

			// compute error estimate:
			mean = value_type(result_kronrod) * 0.5;
			result_asc = wgk[n - 1] * abs(value_type(f_center) - value_type(mean));
			for(j = 0; j < n - 1; ++j){
				result_asc += wgk[j] * (abs(fv1[j]-mean) + abs(fv2[j]-mean));
			}
			err = (result_kronrod - result_gauss) * width;
			result_abs *= width;
			result_asc *= width;
			if(result_asc != 0 && err != 0){
				value_type scale = pow((200 * err / result_asc), 1.5);
				if (scale < 1){
					err = result_asc * scale;
				}else{
					err = result_asc;
				}
			}
			if(result_abs > ScalarTraits<value_type>::min_value() / (50 * ScalarTraits<value_type>::epsilon())){
				value_type min_err = 50 * ScalarTraits<value_type>::epsilon() * result_abs;
				if(min_err > err){
					err = min_err;
				}
			}
			ve.err = err;

			return 0; // no choice but to divide 0th dimension
		}
	};
	int ruleadapt_integrate(Rule *r, const Hypercube &h, return_type &val, value_type *err){
		size_t initialRegions;	// number of initial regions (non-adaptive)
		size_t minIter;		// minimum number of adaptive subdivisions
		size_t maxIter;		// maximum number of adaptive subdivisions
		size_t initialPoints;
		std::priority_queue<Region> regions;
		size_t i;
		int status = -1; // = ERROR

		initialRegions = 1; // or: use some percentage of maxEval/r->num_points
		initialPoints = initialRegions * r->num_points;
		minIter = 0;
		if(_params.max_evals > 0){
			if (initialPoints > _params.max_evals) {
				initialRegions = _params.max_evals / r->num_points;
				initialPoints = initialRegions * r->num_points;
			}
			_params.max_evals -= initialPoints;
			maxIter = _params.max_evals / (2 * r->num_points);
		}else{
			maxIter = (size_t)(-1);
		}

		if(initialRegions == 0){
			return status;	// ERROR
		}

		value_type total_val = 0;
		value_type total_err = 0;
		{
			Region temp_region(h);
			temp_region.Eval(this->f, *r);
			regions.push(temp_region);
			total_val += value_type(temp_region.ve.val);
			total_err += temp_region.ve.err;
		}
		
		int status_update_interval = maxIter/100;
		if(0 == status_update_interval){ status_update_interval = 1; }
		for (i = 0; i < maxIter; ++i) {
#ifdef USE_INTEGRATION_STATUS_CALLBACK
			Cubature::IntegrationStatus Istatus;
#endif
			if (i >= minIter && (total_err <= _params.max_absolute_error || (total_err/total_val) <= _params.max_relative_error)) {
				status = 0; // converged
#ifdef USE_INTEGRATION_STATUS_CALLBACK
				if(NULL != status_func){
					Istatus.num_evals = initialPoints+i*2*r->num_points;
					Istatus.cur_value = total_val;
					Istatus.cur_error = total_err;
					(*status_func)(Istatus);
				}
#endif
				break;
			}
#ifdef USE_INTEGRATION_STATUS_CALLBACK
			if(0 == (i%status_update_interval)){
				if(NULL != status_func){
					Istatus.num_evals = initialPoints+i*2*r->num_points;
					Istatus.cur_value = total_val;
					Istatus.cur_error = total_err;
					(*status_func)(Istatus);
				}
			}
#endif
			Region R = regions.top();
			regions.pop();
			total_val -= value_type(R.ve.val);
			total_err -= R.ve.err;
			Region R2(R);
			R.Split(R2);

			R.Eval(this->f, *r);
			regions.push(R);
			total_val += value_type(R.ve.val);
			total_err += R.ve.err;

			R2.Eval(this->f, *r);
			regions.push(R2);
			total_val += value_type(R2.ve.val);
			total_err += R2.ve.err;
		}

		value_type error;
		while(!regions.empty()){
			val += regions.top().ve.val;
			error += regions.top().ve.err;
			regions.pop();
		}
		if(NULL != err){ *err = error; }

		return status;
	}
private:
	Parameters _params;
	Hypercube h;
};

}; // namespace Integration

#endif // _INTEGRATOR_H_
