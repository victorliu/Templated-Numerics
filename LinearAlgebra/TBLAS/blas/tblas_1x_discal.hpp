template <typename DiagT, typename ScalarT>
void TBLAS_NAME(discal,DISCAL)(
  const size_t &n,
  const DiagT *D, const size_t &incd,
  ScalarT *x, const size_t &incx){
	size_t id = TBLAS_STARTING_OFFSET(n, incd);
	size_t ix = TBLAS_STARTING_OFFSET(n, incx);
	for(size_t i = 0; i < n; ++i){
		x[ix] *= D[id];
		id += incd;
		ix += incx;
	}
}
