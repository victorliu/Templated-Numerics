template <typename ScalarT>
void TBLAS_NAME(rotg,ROTG)(
  ScalarT &a,
  ScalarT &b,
  ScalarT &c,
  ScalarT &s){
	// Taken from zrotg
	if(0 == TBLAS_ABS(a)){
		c = 0;
		s = 1;
		a = b;
		return;
	}
	typename TBLAS_TRAITS(ScalarT)::real_t scale = TBLAS_ABS(a) + TBLAS_ABS(b);
	typename TBLAS_TRAITS(ScalarT)::real_t  as = TBLAS_ABS(a/scale), bs = TBLAS_ABS(b/scale);
	typename TBLAS_TRAITS(ScalarT)::real_t  norm = scale*TBLAS_SQRT(as*as + bs*bs);
	ScalarT alpha = a/TBLAS_ABS(a);
	c = TBLAS_ABS(a)/norm;
	s = alpha * TBLAS_CONJ(b)/norm;
	a = alpha*norm;
}
