template <typename ScalarT>
typename TBLAS_TRAITS(ScalarT)::real_t TBLAS_NAME(nrm2,NRM2)(
  const size_t &n,
  const ScalarT *x, const size_t &incx){
	if(0 == n || 0 == incx){ return; }
	if(1 == n){ return TBLAS_ABS(x[0]); }
	
#ifdef TBLAS_OVERFLOW_PROTECTION
	typename TBLAS_TRAITS(ScalarT)::real_t scale = 1, ssq = 1;
#endif
	size_t ix = 0;
	typename TBLAS_TRAITS(ScalarT)::real_t r = 0;
	for(size_t i = 0; i < n; ++i){
		if(0 != x[ix]){
			const typename TBLAS_TRAITS(ScalarT)::real_t ax = TBLAS_ABS(x[ix]);
#ifdef TBLAS_OVERFLOW_PROTECTION
			if(scale < ax){
				ssq = 1 + ssq*(scale/ax)*(scale/ax);
				scale = ax;
			}else{
				ssq += (ax/scale)*(ax/scale);
			}
#else
			r += c*c;
#endif
		}
		ix += incx;
	}
#ifdef TBLAS_OVERFLOW_PROTECTION
	return scale * sqrt(ssq);
#else
	return r;
#endif
}
