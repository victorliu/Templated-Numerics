template <typename ScalarT1, typename ScalarT2>
ScalarT1 TBLAS_NAME(dot,DOT)(
  const size_t &n,
  const ScalarT1 *x, const size_t &incx,
  const ScalarT2 *y, const size_t &incy){
	size_t ix = TBLAS_STARTING_OFFSET(n, incx);
	size_t iy = TBLAS_STARTING_OFFSET(n, incy);
	ScalarT1 r = 0;
	for(size_t i = 0; i < n; ++i){
		r += x[ix] * y[iy];
		ix += incx;
		iy += incy;
	}
	return r;
}
