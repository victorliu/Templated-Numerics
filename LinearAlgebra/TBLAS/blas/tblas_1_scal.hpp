template <typename ScalarT>
void TBLAS_NAME(scal,SCAL)(
  const size_t &n,
  const ScalarT &alpha,
  ScalarT *x, const size_t &incx){
	if(incx <= 0){ return; }
	size_t ix = TBLAS_STARTING_OFFSET(n, incx);
	for(size_t i = 0; i < n; ++i){
		x[ix] *= alpha;
		ix += incx;
	}
}
