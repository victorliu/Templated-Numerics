template <typename ScalarT>
typename TBLAS_TRAITS(ScalarT)::real_t TBLAS_NAME(iamax,IAMAX)(
  const size_t &n,
  const ScalarT *x, const size_t &incx){
	if(0 == incx){ return 0; }
	
	size_t ix = 0;
	size_t r = 0;
	typename TBLAS_TRAITS(ScalarT)::real_t max = 0;
	for(size_t i = 0; i < n; ++i){
		typename TBLAS_TRAITS(ScalarT)::real_t ax = TBLAS_ABS(x[ix]);
		if(ax > max){
			max = ax;
			r = i;
		}
		ix += incx;
	}
	return r;
}
