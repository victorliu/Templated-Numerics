#ifndef _INTEGRATOR1_H_
#define _INTEGRATOR1_H_

// Preprocessor macros
//   USE_INTEGRATION_STATUS_CALLBACK

#ifdef USE_INTEGRATION_STATUS_CALLBACK
# include "Status.h"
#endif

namespace Integration{

template <class DomainType, class ReturnType>
class Function1{
public:
	typedef DomainType domain_type;
	typedef ReturnType return_type;
	virtual return_type operator()(const domain_type &x) const = 0;
	virtual size_t MaxParallel() const{ return 1; }
};


template <typename FunctionType>
class Integrator1{
public:
	typedef FunctionType function_type;
	typedef typename FunctionType::domain_type domain_type;
	typedef typename FunctionType::return_type return_type;
	
	
	Integrator1(function_type func, const domain_type &a_, const domain_type &b_
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		, const StatusCallback &status_callback = DefaultStatusCallback()
#endif
		):
		f(func), a(a_), b(b_)
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		, status_func(status_callback)
#endif
	{
		if(b < a){ std::swap(a, b); }
	}
	virtual return_type Integrate(domain_type *error = NULL) = 0;
	virtual return_type Refine(domain_type *error = NULL) = 0;
protected:
	FunctionType f;
	domain_type a, b;
#ifdef USE_INTEGRATION_STATUS_CALLBACK
	StatusCallback status_func;
#endif
};

template <typename FunctionType>
class MidpointRule : public Integrator1<FunctionType>{
public:
	typedef FunctionType function_type;
	typedef typename FunctionType::domain_type domain_type;
	typedef typename FunctionType::return_type return_type;
	
	struct Parameters{
		size_t max_evals; // must be at least 2
	};
	MidpointRule(
		FunctionType func, const domain_type &a, const domain_type &b,
		const Parameters &params
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		, const StatusCallback &status_callback = DefaultStatusCallback()
#endif
		):
		Integrator1<FunctionType>(func, a, b
#ifdef USE_INTEGRATION_STATUS_CALLBACK
			, status_callback
#endif
		),
		p(params)
	{}
	return_type Integrate(domain_type *error = NULL){
		domain_type dx = domain_type(this->b-this->a)/domain_type(p.max_evals);
		domain_type hdx = 0.5*dx;
		value = return_type(0);
		domain_type im1 = domain_type(this->b-this->a)/domain_type(p.max_evals-1);
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		size_t status_interval = p.max_evals/status_func.NumUpdates();
#endif
		for(size_t i = 0; i < p.max_evals; ++i){
			value += dx * f(domain_type(i)*im1 + hdx);
#ifdef USE_INTEGRATION_STATUS_CALLBACK
			if(status_func.WantUpdates() && (i > 0) && (0 == (i%status_interval))){
				StatusCallback status;
				status.num_evals = i+1;
				status.current_value = value;
				status.current_error = 0;
				status_func(status);
			}
#endif
		}
		// Estimate error here, if possible
		domain_type err = 0;
		if(NULL != error){ *error = err; }
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		if(status_func.WantUpdates()){
			StatusCallback status;
			status.num_evals = p.max_evals;
			status.current_value = value;
			status.current_error = error;
			status_func(status);
		}
#endif
		return value;
	}
	return_type Refine(domain_type *error = NULL){
		p.max_evals *= 2;
		value *= (return_type(1)/return_type(2));
		domain_type dx = domain_type(this->b-this->a)/domain_type(p.max_evals);
		domain_type hdx = 0.5*dx;
		domain_type im1 = domain_type(this->b-this->a)/domain_type(p.max_evals-1);
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		size_t status_interval = p.max_evals/status_func.NumUpdates();
#endif
		for(size_t i = 0; i < p.max_evals; ++i){
			value += dx * f(domain_type(i)*im1 + hdx);
#ifdef USE_INTEGRATION_STATUS_CALLBACK
			if(status_func.WantUpdates() && (i > 0) && (0 == (i%status_interval))){
				StatusCallback status;
				status.num_evals = i+1;
				status.current_value = value;
				status.current_error = 0;
				status_func(status);
			}
#endif
		}
		// Estimate error here, if possible
		domain_type err = 0;
		if(NULL != error){ *error = err; }
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		if(status_func.WantUpdates()){
			StatusCallback status;
			status.num_evals = p.max_evals;
			status.current_value = value;
			status.current_error = error;
			status_func(status);
		}
#endif
		return value;
	}
private:
	return_type value;
	Parameters p;
};





template <typename FunctionType>
class GaussKronrodRule : public Integrator1<FunctionType>{
public:
	typedef FunctionType function_type;
	typedef typename FunctionType::domain_type domain_type;
	typedef typename FunctionType::return_type return_type;
	
	struct Parameters{
		size_t key;
		size_t max_evals;
		domain_type max_absolute_error;
		domain_type max_relative_error;
	};
	GaussKronrodRule(
		FunctionType func, const domain_type &a, const domain_type &b,
		const Parameters &params
#ifdef USE_INTEGRATION_STATUS_CALLBACK
		, const StatusCallback &status_callback = DefaultStatusCallback()
#endif
		):
		Integrator1<FunctionType>(func, a, b
#ifdef USE_INTEGRATION_STATUS_CALLBACK
			, status_callback
#endif
		),
		p(params)
	{}
	return_type Integrate(domain_type *error = NULL){
		return_type result;
		domain_type abserr;
		size_t neval;
		int ier;
		size_t limit = p.max_evals;
		size_t last;
		int *iwork = new int[limit];
		dqag(this->f, this->a, this->b, p.max_absolute_error, p.max_relative_error, p.key, result, abserr, neval, ier, limit, last, iwork);
		if(NULL != error){
		*error = abserr;
		}
		return result;
	}
	return_type Refine(domain_type *error = NULL){
		return value;
	}
private:
	return_type value;
	Parameters p;
	
	
	
	void dqk15(
		const function_type &f,
		const domain_type &a,
		const domain_type &b,
		// Outputs
		return_type &result,
		domain_type &abserr,
		domain_type &resabs,
		domain_type &resasc
	){
		static const double wg[4] = { .129484966168869693270611432679082,
			.27970539148927666790146777142378,
			.381830050505118944950369775488975,
			.417959183673469387755102040816327 };
		static const double xgk[8] = { .991455371120812639206854697526329,
			.949107912342758524526189684047851,
			.864864423359769072789712788640926,
			.741531185599394439863864773280788,
			.58608723546769113029414483825873,
			.405845151377397166906606412076961,
			.207784955007898467600689403773245,0. };
		static const double wgk[8] = { .02293532201052922496373200805897,
			.063092092629978553290700663189204,
			.104790010322250183839876322541518,
			.140653259715525918745189590510238,
			.16900472663926790282658342659855,
			.190350578064785409913256402421014,
			.204432940075298892414161999234649,
			.209482141084727828012999174891714 };

		return_type fv1[7], fv2[7];

		const domain_type center = 0.5*(a + b);
		const domain_type half_length = 0.5*(b - a);
		const domain_type half_length_abs = std::abs(half_length);

		// compute the 15-point kronrod approximation to
		// the integral, and estimate the absolute error.
		const return_type fc = f(center);
		return_type resg = fc * wg[3];
		return_type resk = fc * wgk[7];
		resabs = std::abs(resk);
		for(size_t j = 0; j < 3; ++j){
			size_t jtw = 2*j + 1;
			domain_type absc = half_length * xgk[jtw];
			fv1[jtw] = f(center - absc);
			fv2[jtw] = f(center + absc);
			return_type fsum = fv1[jtw] + fv2[jtw];
			resg += wg[j] * fsum;
			resk += wgk[jtw] * fsum;
			resabs += wgk[jtw] * (std::abs(fv1[jtw]) + std::abs(fv2[jtw]));
		}
		for(size_t j = 0; j < 4; ++j){
			size_t jtwm1 = 2*j;
			domain_type absc = half_length * xgk[jtwm1];
			fv1[jtwm1] = f(center - absc);
			fv2[jtwm1] = f(center + absc);
			return_type fsum = fv1[jtwm1] + fv2[jtwm1];
			resk += wgk[jtwm1] * fsum;
			resabs += wgk[jtwm1] * (std::abs(fv1[jtwm1]) + std::abs(fv2[jtwm1]));
		}
		return_type reskh = resk * .5;
		resasc = wgk[7] * std::abs(fc - reskh);
		for(size_t j = 0; j < 7; ++j){
			resasc += wgk[j] * (std::abs(fv1[j] - reskh) + std::abs(fv2[j] - reskh));
		}
		result = resk * half_length;
		resabs *= half_length_abs;
		resasc *= half_length_abs;
		abserr = std::abs((resk - resg) * half_length);
		if(resasc != 0. && abserr != 0.){
			abserr = resasc * std::min(1.0,pow(abserr * 200. / resasc, 1.5));
		}
		const domain_type epmach = 2*std::numeric_limits<domain_type>::epsilon();
		if (resabs > std::numeric_limits<domain_type>::min() / (epmach * 50.)){
			abserr = std::max(epmach * 50. * resabs, abserr);
		}
		return;
	}
	
	void dqk21(
		const function_type &f,
		const domain_type &a,
		const domain_type &b,
		// Outputs
		return_type &result,
		domain_type &abserr,
		domain_type &resabs,
		domain_type &resasc
	){
		static double wg[5] = { .066671344308688137593568809893332,
			.149451349150580593145776339657697,
			.219086362515982043995534934228163,
			.269266719309996355091226921569469,
			.295524224714752870173892994651338 };
		static double xgk[11] = { .995657163025808080735527280689003,
			.973906528517171720077964012084452,
			.930157491355708226001207180059508,
			.865063366688984510732096688423493,
			.780817726586416897063717578345042,
			.679409568299024406234327365114874,
			.562757134668604683339000099272694,
			.433395394129247190799265943165784,
			.294392862701460198131126603103866,
			.14887433898163121088482600112972,0. };
		static double wgk[11] = { .011694638867371874278064396062192,
			.03255816230796472747881897245939,
			.05475589657435199603138130024458,
			.07503967481091995276704314091619,
			.093125454583697605535065465083366,
			.109387158802297641899210590325805,
			.123491976262065851077958109831074,
			.134709217311473325928054001771707,
			.142775938577060080797094273138717,
			.147739104901338491374841515972068,
			.149445554002916905664936468389821 };


		return_type fv1[10], fv2[10];
		return_type resg, resk, fsum, fval1, fval2;

		const domain_type epmach = 2*std::numeric_limits<domain_type>::epsilon();
		const domain_type uflow = std::numeric_limits<domain_type>::min();

		const domain_type center = 0.5*(a + b);
		const domain_type half_length = 0.5*(b - a);
		const domain_type half_length_abs = std::abs(half_length);

		// compute the 21-point kronrod approximation to
		// the integral, and estimate the absolute error.

		resg = 0.;
		const return_type fc = f(center);
		resk = wgk[10] * fc;
		resabs = std::abs(resk);
		for(size_t j = 0; j < 5; ++j) {
			size_t jtw = 2*j+1;
			domain_type absc = half_length * xgk[jtw];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtw] = fval1;
			fv2[jtw] = fval2;
			fsum = fval1 + fval2;
			resg += wg[j] * fsum;
			resk += wgk[jtw] * fsum;
			resabs += wgk[jtw] * (std::abs(fval1) + std::abs(fval2));
		}
		for(size_t j = 0; j < 5; ++j) {
			size_t jtwm1 = 2*j;
			domain_type absc = half_length * xgk[jtwm1];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtwm1] = fval1;
			fv2[jtwm1] = fval2;
			fsum = fval1 + fval2;
			resk += wgk[jtwm1] * fsum;
			resabs += wgk[jtwm1] * (std::abs(fval1) + std::abs(fval2));
		}
		return_type reskh = resk * .5;
		resasc = wgk[10] * std::abs(fc - reskh);
		for(size_t j = 0; j < 10; ++j) {
			resasc += wgk[j] * (std::abs(fv1[j] - reskh) + (std::abs(fv2[j] - reskh)));
		}
		result = resk * half_length;
		resabs *= half_length_abs;
		resasc *= half_length_abs;
		abserr = std::abs((resk - resg) * half_length);
		if (resasc != 0. && abserr != 0.) {
			abserr = resasc * std::min(1.0, pow(abserr * 200. / resasc, 1.5));
		}
		if (resabs > uflow / (epmach * 50.)) {
			abserr = std::max(epmach * 50. * resabs, abserr);
		}
	}

	void dqk31(
		const function_type &f,
		const domain_type &a,
		const domain_type &b,
		// Outputs
		return_type &result,
		domain_type &abserr,
		domain_type &resabs,
		domain_type &resasc
	){
		static const double wg[8] = { .030753241996117268354628393577204,
			.070366047488108124709267416450667,
			.107159220467171935011869546685869,
			.139570677926154314447804794511028,
			.166269205816993933553200860481209,
			.186161000015562211026800561866423,
			.198431485327111576456118326443839,
			.202578241925561272880620199967519 };
		static const double xgk[16] = { .998002298693397060285172840152271,
			.987992518020485428489565718586613,
			.967739075679139134257347978784337,
			.937273392400705904307758947710209,
			.897264532344081900882509656454496,
			.848206583410427216200648320774217,
			.790418501442465932967649294817947,
			.724417731360170047416186054613938,
			.650996741297416970533735895313275,
			.570972172608538847537226737253911,
			.485081863640239680693655740232351,
			.394151347077563369897207370981045,
			.299180007153168812166780024266389,
			.201194093997434522300628303394596,
			.101142066918717499027074231447392,0. };
		static const double wgk[16] = { .005377479872923348987792051430128,
			.015007947329316122538374763075807,
			.025460847326715320186874001019653,
			.03534636079137584622203794847836,
			.04458975132476487660822729937328,
			.05348152469092808726534314723943,
			.062009567800670640285139230960803,
			.069854121318728258709520077099147,
			.076849680757720378894432777482659,
			.083080502823133021038289247286104,
			.088564443056211770647275443693774,
			.093126598170825321225486872747346,
			.096642726983623678505179907627589,
			.099173598721791959332393173484603,
			.10076984552387559504494666261757,
			.101330007014791549017374792767493 };

		return_type fv1[15], fv2[15];
		return_type resg, resk, fsum, fval1, fval2;

		const domain_type epmach = 2*std::numeric_limits<domain_type>::epsilon();
		const domain_type uflow = std::numeric_limits<domain_type>::min();

		const domain_type center = 0.5*(a + b);
		const domain_type half_length = 0.5*(b - a);
		const domain_type half_length_abs = std::abs(half_length);

	// compute the 31-point kronrod approximation to
	// the integral, and estimate the absolute error.

		const return_type fc = f(center);
		resg = wg[7] * fc;
		resk = wgk[15] * fc;
		resabs = std::abs(resk);
		for(size_t j = 0; j < 7; ++j) {
			size_t jtw = 2*j+1;
			domain_type absc = half_length * xgk[jtw];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtw] = fval1;
			fv2[jtw] = fval2;
			fsum = fval1 + fval2;
			resg += wg[j] * fsum;
			resk += wgk[jtw] * fsum;
			resabs += wgk[jtw] * (std::abs(fval1) + std::abs(fval2));
		}
		for(size_t j = 0; j < 8; ++j) {
			size_t jtwm1 = 2*j;
			domain_type absc = half_length * xgk[jtwm1];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtwm1] = fval1;
			fv2[jtwm1] = fval2;
			fsum = fval1 + fval2;
			resk += wgk[jtwm1] * fsum;
			resabs += wgk[jtwm1] * (std::abs(fval1) + std::abs(fval2));
		}
		return_type reskh = resk * .5;
		resasc = wgk[15] * std::abs(fc - reskh);
		for (size_t j = 0; j < 15; ++j) {
			resasc += wgk[j] * (std::abs(fv1[j] - reskh) + std::abs(fv2[j] - reskh));
		}
		result = resk * half_length;
		resabs *= half_length_abs;
		resasc *= half_length_abs;
		abserr = std::abs((resk - resg) * half_length);
		if (resasc != 0. && abserr != 0.) {
			abserr = resasc * std::min(1.0,pow(abserr * 200. / resasc, 1.5));
		}
		if (resabs > uflow / (epmach * 50.)) {
			abserr = std::max(epmach * 50. * resabs,abserr);
		}
	}

	void dqk41(
		const function_type &f,
		const domain_type &a,
		const domain_type &b,
		// Outputs
		return_type &result,
		domain_type &abserr,
		domain_type &resabs,
		domain_type &resasc
	){
		static double wg[10] = { .017614007139152118311861962351853,
			.040601429800386941331039952274932,
			.062672048334109063569506535187042,
			.083276741576704748724758143222046,
			.10193011981724043503675013548035,
			.118194531961518417312377377711382,
			.131688638449176626898494499748163,
			.142096109318382051329298325067165,
			.149172986472603746787828737001969,
			.152753387130725850698084331955098 };
		static double xgk[21] = { .998859031588277663838315576545863,
			.99312859918509492478612238847132,
			.981507877450250259193342994720217,
			.963971927277913791267666131197277,
			.940822633831754753519982722212443,
			.912234428251325905867752441203298,
			.878276811252281976077442995113078,
			.839116971822218823394529061701521,
			.795041428837551198350638833272788,
			.746331906460150792614305070355642,
			.693237656334751384805490711845932,
			.636053680726515025452836696226286,
			.575140446819710315342946036586425,
			.510867001950827098004364050955251,
			.44359317523872510319999221349264,
			.373706088715419560672548177024927,
			.301627868114913004320555356858592,
			.227785851141645078080496195368575,
			.152605465240922675505220241022678,
			.076526521133497333754640409398838,0. };
		static double wgk[21] = { .003073583718520531501218293246031,
			.008600269855642942198661787950102,
			.014626169256971252983787960308868,
			.020388373461266523598010231432755,
			.025882133604951158834505067096153,
			.031287306777032798958543119323801,
			.036600169758200798030557240707211,
			.041668873327973686263788305936895,
			.046434821867497674720231880926108,
			.050944573923728691932707670050345,
			.055195105348285994744832372419777,
			.059111400880639572374967220648594,
			.062653237554781168025870122174255,
			.065834597133618422111563556969398,
			.068648672928521619345623411885368,
			.07105442355344406830579036172321,
			.073030690332786667495189417658913,
			.074582875400499188986581418362488,
			.075704497684556674659542775376617,
			.076377867672080736705502835038061,
			.076600711917999656445049901530102 };


		return_type fv1[20], fv2[20];
		return_type resg, resk, fsum, fval1, fval2;

		const domain_type epmach = 2*std::numeric_limits<domain_type>::epsilon();
		const domain_type uflow = std::numeric_limits<domain_type>::min();

		const domain_type center = 0.5*(a + b);
		const domain_type half_length = 0.5*(b - a);
		const domain_type half_length_abs = std::abs(half_length);

		// compute the 41-point gauss-kronrod approximation to
		// the integral, and estimate the absolute error.

		resg = 0.;
		const return_type fc = f(center);
		resk = wgk[20] * fc;
		resabs = std::abs(resk);
		for(size_t j = 0; j < 10; ++j) {
			size_t jtw = 2*j+1;
			domain_type absc = half_length * xgk[jtw - 1];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtw] = fval1;
			fv2[jtw] = fval2;
			fsum = fval1 + fval2;
			resg += wg[j] * fsum;
			resk += wgk[jtw] * fsum;
			resabs += wgk[jtw] * (std::abs(fval1) + std::abs(fval2));
		}
		for(size_t j = 0; j < 10; ++j) {
			size_t jtwm1 = 2*j;
			domain_type absc = half_length * xgk[jtwm1];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtwm1] = fval1;
			fv2[jtwm1] = fval2;
			fsum = fval1 + fval2;
			resk += wgk[jtwm1] * fsum;
			resabs += wgk[jtwm1] * (std::abs(fval1) + std::abs(fval2));
		}
		return_type reskh = resk * .5;
		resasc = wgk[20] * std::abs(fc - reskh);
		for(size_t j = 1; j <= 20; ++j) {
			resasc += wgk[j - 1] * (std::abs(fv1[j - 1] - reskh) + std::abs(fv2[j - 1] - reskh));
		}
		result = resk * half_length;
		resabs *= half_length_abs;
		resasc *= half_length_abs;
		abserr = std::abs((resk - resg) * half_length);
		if (resasc != 0. && abserr != 0.) {
			abserr = resasc * std::min(1.0, pow(abserr * 200. / resasc, 1.5));
		}
		if (resabs > uflow / (epmach * 50.)) {
			abserr = std::max(epmach * 50. * resabs, abserr);
		}
	}


	void dqk51(
		const function_type &f,
		const domain_type &a,
		const domain_type &b,
		// Outputs
		return_type &result,
		domain_type &abserr,
		domain_type &resabs,
		domain_type &resasc
	){
		static double wg[13] = { .011393798501026287947902964113235,
			.026354986615032137261901815295299,
			.040939156701306312655623487711646,
			.054904695975835191925936891540473,
			.068038333812356917207187185656708,
			.080140700335001018013234959669111,
			.091028261982963649811497220702892,
			.100535949067050644202206890392686,
			.108519624474263653116093957050117,
			.114858259145711648339325545869556,
			.119455763535784772228178126512901,
			.122242442990310041688959518945852,
			.12317605372671545120390287307905 };
		static double xgk[26] = { .999262104992609834193457486540341,
			.995556969790498097908784946893902,
			.988035794534077247637331014577406,
			.976663921459517511498315386479594,
			.961614986425842512418130033660167,
			.942974571228974339414011169658471,
			.920747115281701561746346084546331,
			.894991997878275368851042006782805,
			.86584706529327559544899696958834,
			.83344262876083400142102110869357,
			.797873797998500059410410904994307,
			.759259263037357630577282865204361,
			.717766406813084388186654079773298,
			.673566368473468364485120633247622,
			.626810099010317412788122681624518,
			.577662930241222967723689841612654,
			.52632528433471918259962377815801,
			.473002731445714960522182115009192,
			.417885382193037748851814394594572,
			.361172305809387837735821730127641,
			.303089538931107830167478909980339,
			.243866883720988432045190362797452,
			.183718939421048892015969888759528,
			.122864692610710396387359818808037,
			.061544483005685078886546392366797,0. };
		static double wgk[26] = { .001987383892330315926507851882843,
			.005561932135356713758040236901066,
			.009473973386174151607207710523655,
			.013236229195571674813656405846976,
			.016847817709128298231516667536336,
			.020435371145882835456568292235939,
			.024009945606953216220092489164881,
			.027475317587851737802948455517811,
			.030792300167387488891109020215229,
			.034002130274329337836748795229551,
			.03711627148341554356033062536762,
			.040083825504032382074839284467076,
			.042872845020170049476895792439495,
			.04550291304992178890987058475266,
			.047982537138836713906392255756915,
			.05027767908071567196332525943344,
			.052362885806407475864366712137873,
			.054251129888545490144543370459876,
			.055950811220412317308240686382747,
			.057437116361567832853582693939506,
			.058689680022394207961974175856788,
			.059720340324174059979099291932562,
			.060539455376045862945360267517565,
			.061128509717053048305859030416293,
			.061471189871425316661544131965264,
			.061580818067832935078759824240066 };


		return_type fv1[25], fv2[25];
		return_type resg, resk, fsum, fval1, fval2;

		const domain_type epmach = 2*std::numeric_limits<domain_type>::epsilon();
		const domain_type uflow = std::numeric_limits<domain_type>::min();

		const domain_type center = 0.5*(a + b);
		const domain_type half_length = 0.5*(b - a);
		const domain_type half_length_abs = std::abs(half_length);

		// compute the 51-point kronrod approximation to
		// the integral, and estimate the absolute error.

		const return_type fc = f(center);
		resg = wg[12] * fc;
		resk = wgk[25] * fc;
		resabs = std::abs(resk);
		for(size_t j = 0; j < 12; ++j) {
			size_t jtw = 2*j+1;
			domain_type absc = half_length * xgk[jtw];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtw] = fval1;
			fv2[jtw] = fval2;
			fsum = fval1 + fval2;
			resg += wg[j] * fsum;
			resk += wgk[jtw] * fsum;
			resabs += wgk[jtw] * (std::abs(fval1) + std::abs(fval2));
		}
		for(size_t j = 0; j < 13; ++j) {
			size_t jtwm1 = 2*j;
			domain_type absc = half_length * xgk[jtwm1];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtwm1] = fval1;
			fv2[jtwm1] = fval2;
			fsum = fval1 + fval2;
			resk += wgk[jtwm1] * fsum;
			resabs += wgk[jtwm1] * (std::abs(fval1) + std::abs(fval2));
		}
		return_type reskh = resk * .5;
		resasc = wgk[25] * std::abs(fc - reskh);
		for(size_t j = 0; j < 25; ++j) {
			resasc += wgk[j] * (std::abs(fv1[j] - reskh) + std::abs(fv2[j] - reskh));
		}
		result = resk * half_length;
		resabs *= half_length_abs;
		resasc *= half_length_abs;
		abserr = std::abs((resk - resg) * half_length);
		if (resasc != 0. && abserr != 0.) {
			abserr = resasc * std::min(1.0, pow(abserr * 200. / resasc, 1.5));
		}
		if (resabs > uflow / (epmach * 50.)) {
			abserr = std::max(epmach * 50. * resabs, abserr);
		}
	}


	void dqk61(
		const function_type &f,
		const domain_type &a,
		const domain_type &b,
		// Outputs
		return_type &result,
		domain_type &abserr,
		domain_type &resabs,
		domain_type &resasc
	){
		static double wg[15] = { .007968192496166605615465883474674,
			.018466468311090959142302131912047,
			.028784707883323369349719179611292,
			.038799192569627049596801936446348,
			.048402672830594052902938140422808,
			.057493156217619066481721689402056,
			.065974229882180495128128515115962,
			.073755974737705206268243850022191,
			.08075589522942021535469493846053,
			.086899787201082979802387530715126,
			.092122522237786128717632707087619,
			.09636873717464425963946862635181,
			.099593420586795267062780282103569,
			.101762389748405504596428952168554,
			.102852652893558840341285636705415 };
		static double xgk[31] = { .999484410050490637571325895705811,
			.996893484074649540271630050918695,
			.991630996870404594858628366109486,
			.983668123279747209970032581605663,
			.973116322501126268374693868423707,
			.960021864968307512216871025581798,
			.944374444748559979415831324037439,
			.926200047429274325879324277080474,
			.905573307699907798546522558925958,
			.882560535792052681543116462530226,
			.857205233546061098958658510658944,
			.829565762382768397442898119732502,
			.799727835821839083013668942322683,
			.767777432104826194917977340974503,
			.733790062453226804726171131369528,
			.69785049479331579693229238802664,
			.660061064126626961370053668149271,
			.620526182989242861140477556431189,
			.57934523582636169175602493217254,
			.536624148142019899264169793311073,
			.492480467861778574993693061207709,
			.447033769538089176780609900322854,
			.400401254830394392535476211542661,
			.352704725530878113471037207089374,
			.304073202273625077372677107199257,
			.254636926167889846439805129817805,
			.204525116682309891438957671002025,
			.153869913608583546963794672743256,
			.102806937966737030147096751318001,
			.051471842555317695833025213166723,0. };
		static double wgk[31] = { .00138901369867700762455159122676,
			.003890461127099884051267201844516,
			.00663070391593129217331982636975,
			.009273279659517763428441146892024,
			.011823015253496341742232898853251,
			.01436972950704580481245143244358,
			.016920889189053272627572289420322,
			.019414141193942381173408951050128,
			.021828035821609192297167485738339,
			.024191162078080601365686370725232,
			.026509954882333101610601709335075,
			.028754048765041292843978785354334,
			.030907257562387762472884252943092,
			.032981447057483726031814191016854,
			.034979338028060024137499670731468,
			.036882364651821229223911065617136,
			.038678945624727592950348651532281,
			.040374538951535959111995279752468,
			.04196981021516424614714754128597,
			.043452539701356069316831728117073,
			.044814800133162663192355551616723,
			.046059238271006988116271735559374,
			.047185546569299153945261478181099,
			.048185861757087129140779492298305,
			.049055434555029778887528165367238,
			.049795683427074206357811569379942,
			.050405921402782346840893085653585,
			.050881795898749606492297473049805,
			.051221547849258772170656282604944,
			.051426128537459025933862879215781,
			.051494729429451567558340433647099 };


		return_type fv1[30], fv2[30];
		return_type resg, resk, fsum, fval1, fval2;

		const domain_type epmach = 2*std::numeric_limits<domain_type>::epsilon();
		const domain_type uflow = std::numeric_limits<domain_type>::min();

		const domain_type center = 0.5*(a + b);
		const domain_type half_length = 0.5*(b - a);
		const domain_type half_length_abs = std::abs(half_length);

		// compute the 61-point kronrod approximation to the
		// integral, and estimate the absolute error.

		resg = 0.;
		const return_type fc = f(center);
		resk = wgk[30] * fc;
		resabs = std::abs(resk);
		for(size_t j = 0; j < 15; ++j) {
			size_t jtw = 2*j+1;
			domain_type absc = half_length * xgk[jtw];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtw] = fval1;
			fv2[jtw] = fval2;
			fsum = fval1 + fval2;
			resg += wg[j] * fsum;
			resk += wgk[jtw] * fsum;
			resabs += wgk[jtw] * (std::abs(fval1) + std::abs(fval2));
		}
		for(size_t j = 0; j < 15; ++j) {
			size_t jtwm1 = 2*j;
			domain_type absc = half_length * xgk[jtwm1];
			fval1 = f(center - absc);
			fval2 = f(center + absc);
			fv1[jtwm1] = fval1;
			fv2[jtwm1] = fval2;
			fsum = fval1 + fval2;
			resk += wgk[jtwm1] * fsum;
			resabs += wgk[jtwm1] * (std::abs(fval1) + std::abs(fval2));
		}
		return_type reskh = resk * .5;
		resasc = wgk[30] * std::abs(fc - reskh);
		for(size_t j = 0; j < 30; ++j) {
			resasc += wgk[j] * (std::abs(fv1[j] - reskh) + std::abs(fv2[j] - reskh));
		}
		result = resk * half_length;
		resabs *= half_length_abs;
		resasc *= half_length_abs;
		abserr = std::abs((resk - resg) * half_length);
		if (resasc != 0. && abserr != 0.) {
			abserr = resasc * std::min(1.0, pow(abserr * 200. / resasc, 1.5));
		}
		if (resabs > uflow / (epmach * 50.)) {
			abserr = std::max(epmach * 50. * resabs, abserr);
		}
	}




	// input: maxerr is the index of the (possibly new) largest error in elist
	//        elist[last] is the index of the new smallest error in elist
	//        iord[last] is not yet defined
	// output: iord contains indices such that elist[iord[0]] ... elist[iord[k]]
	//         is a decreasing sequence with k=last if last<=limit/2+1, otherwise k=limit+1-last
	// Basically, calling this function will insert two error estimates into elist;
	// The larger of the two is contained in elist[maxerr], the smaller is at elist[last]
	void dqpsrt(
		int limit, // maximum number of elements that elist can contain (allocation size)
		int last, // current number of elements in elist and iord
		int &maxerr, // index of nrmax-th largest error, maxerr=iord[nrmax]
		domain_type &ermax, // ermax = elist[maxerr] = elist[iord[nrmax]]
		domain_type *elist, // list of error estimates
		int *iord // list of indices maintaining the sort order
	){
		if(last < 2){
			iord[0] = 0;
			iord[1] = 1;
		}else{
			// this part of the routine is only executed if, due to a
			// difficult integrand, subdivision increased the error
			// estimate. in the normal case the insert procedure should
			// start after the nrmax-th largest error estimate.
			domain_type errmax = elist[maxerr];

			// compute the number of elements in the list to be maintained
			// in descending order. this number depends on the number of
			// subdivisions still allowed.
			size_t jupbn = last; // jupbn is upper bound for j; highest index we can use in elist and iord
			if(last > limit / 2 + 1){
				jupbn = limit - 1 - last;
			}
			domain_type errmin = elist[last]; // the newest smaller error estimate

			// insert errmax by traversing the list top-down,
			// starting comparison from the element elist(iord(nrmax+1)).
			size_t jbnd = jupbn - 1; // avoid the last slot
			bool found_max = false;
			size_t i;
			if(jbnd > 0){ // 0 is where maxerr is right now, so no need to check it
				for(i = 1; i <= jbnd; ++i){
					// Bubble maxerr entry starting from the largest side
					size_t isucc = iord[i];
					if(errmax >= elist[isucc]){
						iord[i-1] = maxerr; // if the position for maxerr was found above, then insert it here
						found_max = true;
						break;
					}
					iord[i-1] = isucc;
				}
			}
			if(!found_max){
				iord[jbnd] = maxerr; // we would get here if maxerr is in fact smaller than all existing entries in elist (except for the last slot). Here we insert it into the last slot
				iord[jupbn] = last; // so the two new error estimates are actually the two smallest ones
			}else{
				// insert errmin by traversing the list bottom-up.
				bool found_min = false;
				for(size_t k = jbnd; k >= i; --k){ // here i is the position right after where maxerr was inserted
					size_t isucc = iord[k];
					if(errmin < elist[isucc]){
						iord[k+1] = last;
						found_min = true;
						break;
					}
					iord[k+1] = isucc;
				}
				if(!found_min){
					iord[i] = last;
				}
			}
		}
		maxerr = iord[0];
		ermax = elist[maxerr];
	}


	void dqage(
		const function_type &f,
		const domain_type &a,
		const domain_type &b,
		const domain_type &epsabs,
		const domain_type &epsrel,
		int key,
		size_t limit,
		// Outputs
		return_type &result,
		domain_type &abserr,
		size_t &neval,
		int &ier,
		domain_type *alist, domain_type *blist, return_type *rlist,
		domain_type *elist, int *iord,
		size_t &last
	){
		domain_type a1, a2, b1, b2;
		return_type area;

		return_type area1, area2, area12;
		domain_type erro12, defab1, defab2;
		int iroff1, iroff2;
		domain_type error1, error2, defabs;
		domain_type errbnd, resabs, errmax;
		int maxerr;
		domain_type errsum;

		const domain_type epmach = 2*std::numeric_limits<domain_type>::epsilon();
		const domain_type uflow = std::numeric_limits<domain_type>::min();

		// test on validity of parameters

		ier = 0;
		neval = 0;
		last = 0;
		result = 0.;
		abserr = 0.;
		alist[0] = a;
		blist[0] = b;
		rlist[0] = 0.;
		elist[0] = 0.;
		iord[0] = 0;
		// Computing MAX
		if(epsabs <= 0. && epsrel < std::max(epmach * 50, 5e-29)){
			ier = 6;
		}
		if(ier == 6){
			return;
		}

		// first approximation to the integral
		if(key <= 0){
			key = 1;
		}
		if(key >= 7){
			key = 6;
		}
		neval = 0;
		switch(key){
		case 1: dqk15(f, a, b, result, abserr, defabs, resabs); break;
		case 2: dqk21(f, a, b, result, abserr, defabs, resabs); break;
		case 3: dqk31(f, a, b, result, abserr, defabs, resabs); break;
		case 4: dqk41(f, a, b, result, abserr, defabs, resabs); break;
		case 5: dqk51(f, a, b, result, abserr, defabs, resabs); break;
		case 6: dqk61(f, a, b, result, abserr, defabs, resabs); break;
		default: break;
		}
		last = 0;
		rlist[0] = result;
		elist[0] = abserr;
		iord[0] = 0;

		// test on accuracy.
		errbnd = std::max(epsabs,epsrel * std::abs(result));
		if (abserr <= epmach * 50. * defabs && abserr > errbnd) {
			ier = 2;
		}
		if (limit == 1) {
			ier = 1;
		}
		if (ier != 0 || (abserr <= errbnd && abserr != resabs) || abserr == 0.){
			if (key != 1) {
				neval = (key * 10 + 1) * (2*neval + 1);
			}
			if (key == 1) {
				neval = neval * 30 + 15;
			}
			return;
		}
		// initialization
		errmax = abserr;
		maxerr = 0;
		area = result;
		errsum = abserr;
		iroff1 = 0;
		iroff2 = 0;

		for(last = 1; last < limit; ++last){
			// bisect the subinterval with the largest error estimate.
			a1 = alist[maxerr];
			b1 = (alist[maxerr] + blist[maxerr]) * .5;
			a2 = b1;
			b2 = blist[maxerr];

			switch(key){
			case 1:
				dqk15(f, a1, b1, area1, error1, resabs, defab1);
				dqk15(f, a2, b2, area2, error2, resabs, defab2);
				break;
			case 2:
				dqk21(f, a1, b1, area1, error1, resabs, defab1);
				dqk21(f, a2, b2, area2, error2, resabs, defab2);
				break;
			case 3:
				dqk31(f, a1, b1, area1, error1, resabs, defab1);
				dqk31(f, a2, b2, area2, error2, resabs, defab2);
				break;
			case 4:
				dqk41(f, a1, b1, area1, error1, resabs, defab1);
				dqk41(f, a2, b2, area2, error2, resabs, defab2);
				break;
			case 5:
				dqk51(f, a1, b1, area1, error1, resabs, defab1);
				dqk51(f, a2, b2, area2, error2, resabs, defab2);
				break;
			case 6:
				dqk61(f, a1, b1, area1, error1, resabs, defab1);
				dqk61(f, a2, b2, area2, error2, resabs, defab2);
				break;
			default: break;
			}

			// improve previous approximations to integral and error and test for accuracy.
			++neval;
			area12 = area1 + area2;
			erro12 = error1 + error2;
			errsum += erro12 - errmax;
			area = area + area12 - rlist[maxerr];
			if(!(defab1 == error1 || defab2 == error2)){
				if (std::abs(rlist[maxerr] - area12) <= std::abs(area12) * 1e-5 && erro12 >= errmax * .99) {
					++iroff1;
				}
				if (last > 9 && erro12 > errmax) {
					++iroff2;
				}
			}
			rlist[maxerr] = area1;
			rlist[last] = area2;

			errbnd = std::max(epsabs,epsrel * std::abs(area));
			if(errsum > errbnd) {
				// test for roundoff error and eventually set error flag.
				if (iroff1 >= 6 || iroff2 >= 20) {
					ier = 2;
				}
				// set error flag in the case that the number of subintervals equals limit.
				if (last == limit) {
					ier = 1;
				}
				// set error flag in the case of bad integrand behaviour at a point of the integration range.
				if(std::max(std::abs(a1),std::abs(b2)) <= (epmach * 100. + 1.) * (std::abs(a2) + uflow * 1e3)){
					ier = 3;
				}
			}
			// append the newly-created intervals to the list.
			if (error2 > error1){
				alist[maxerr] = a2;
				alist[last] = a1;
				blist[last] = b1;
				rlist[maxerr] = area2;
				rlist[last] = area1;
				elist[maxerr] = error2;
				elist[last] = error1;
			}else{
				alist[last] = a2;
				blist[maxerr] = b1;
				blist[last] = b2;
				elist[maxerr] = error1;
				elist[last] = error2;
			}

			// call subroutine dqpsrt to maintain the descending ordering
			// in the list of error estimates and select the subinterval
			// with the largest error estimate (to be bisected next).
			dqpsrt(limit, last, maxerr, errmax, elist, iord);

			// jump out of do-loop
			if (ier != 0 || errsum <= errbnd) {
				break;
			}
		}

		// compute final result.
		result = 0.;
		for(size_t k = 0; k <= last; ++k) {
			result += rlist[k];
		}
		abserr = errsum;

		if (key != 1) {
			neval = (key * 10 + 1) * (2*neval + 1);
		}
		if (key == 1) {
			neval = neval * 30 + 15;
		}
		return;
	}


	void dqag(
		const function_type &f,
		const domain_type &a,
		const domain_type &b,
		const domain_type &epsabs,
		const domain_type &epsrel,
		int key,
		// outputs
		return_type &result, 
		domain_type &abserr,
		size_t &neval,
		int &ier,
		// Dimensions
		size_t limit,
		size_t &last,
		// Work arrays
		int *iwork
	){
		ier = 6;
		neval = 0;
		last = 0;
		result = 0.;
		abserr = 0.;
		if(limit < 1){
			return;
		}
		domain_type *alist = new domain_type[limit];
		domain_type *blist = new domain_type[limit];
		return_type *rlist = new return_type[limit];
		domain_type *elist = new domain_type[limit];
		dqage(f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iwork, last);
		delete [] alist;
		delete [] blist;
		delete [] rlist;
		delete [] elist;
	}

};

}; // namespace Integration

#endif // _INTEGRATOR1_H_
