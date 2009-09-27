template <typename DiagT, typename ScalarT>
void TBLAS_NAME(discal,DISCAL)(
  const size_t &n,
  const DiagT *D, const size_t &incd,
  const ScalarT *x, const size_t &incx,
  const ScalarT &beta,
  ScalarT *y, const size_t &incy){
	size_t id = TBLAS_STARTING_OFFSET(n, incd);
	size_t ix = TBLAS_STARTING_OFFSET(n, incx);
	size_t iy = TBLAS_STARTING_OFFSET(n, incy);
	for(size_t i = 0; i < n; ++i){
		y[iy] = D[id] * x[ix] + beta * y[iy];
		id += incd;
		ix += incx;
		iy += incy;
	}
}
