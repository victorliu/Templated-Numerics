template <typename ScalarT>
void TBLAS_NAME(axpy,AXPY)(
  const size_t &n,
  const ScalarT &alpha,
  const ScalarT *x, const size_t &incx,
  ScalarT *y, const size_t &incy){
	size_t ix = TBLAS_STARTING_OFFSET(n, incx);
	size_t iy = TBLAS_STARTING_OFFSET(n, incy);
	for(size_t i = 0; i < n; ++i){
		y[iy] += alpha * x[ix];
		ix += incx;
		iy += incy;
	}
}
