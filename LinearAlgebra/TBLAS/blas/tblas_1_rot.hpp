template <typename ScalarRealT>
void TBLAS_NAME(rot,ROT)(
  const size_t &n,
  ScalarRealT *x, const size_t &incx,
  ScalarRealT *y, const size_t &incy,
  const ScalarRealT &c, const ScalarRealT &s){
	size_t ix = TBLAS_STARTING_OFFSET(n, incx);
	size_t iy = TBLAS_STARTING_OFFSET(n, incy);
	for(size_t i = 0; i < n; ++i){
		const ScalarRealT xx = x[ix];
		const ScalarRealT yy = x[iy];
		x[ix] = c*xx + s*yy;
		y[iy] = -s*xx + c*yy;
		ix += incx;
		iy += incy;
	}
}
