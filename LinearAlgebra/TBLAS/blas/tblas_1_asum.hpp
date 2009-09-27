template <typename ScalarT>
typename TBLAS_TRAITS(ScalarT)::real_t TBLAS_NAME(asum,ASUM)(
  const size_t &n,
  const ScalarT *x, const size_t &incx){
	if(0 == n || 0 == incx){ return 0; }
	if(1 == n){ return TBLAS_ABS(x[0]); }
	
	size_t ix = 0;
	typename TBLAS_TRAITS(ScalarT)::real_t r = 0;
	for(size_t i = 0; i < n; ++i){
		if(0 != x[ix]){
			r += TBLAS_ABS(x[ix]);
		}
		ix += incx;
	}
	return r;
}
