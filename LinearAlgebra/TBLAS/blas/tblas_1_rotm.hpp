template <typename ScalarRealT>
void TBLAS_NAME(rotm,ROTM)(
  const size_t &n,
  ScalarRealT *x, const size_t &incx,
  ScalarRealT *y, const size_t &incy,
  const ScalarRealT param[5]){
	ScalarRealT h11, h21, h12, h22;

	if(param[0] == -1){
		h11 = param[1];
		h21 = param[2];
		h12 = param[3];
		h22 = param[4];
	}else if(param[0] == 0){
		h11 = 1.0;
		h21 = param[2];
		h12 = param[3];
		h22 = 1.0;
	}else if(param[0] == 1){
		h11 = param[1];
		h21 = -1.0;
		h12 = 1.0;
		h22 = param[4];
	}else if(param[0] == -2){
		return;
	}else{
		TBLAS_ERROR("unrecognized value of param[0]");
		return;
	}

	size_t ix = TBLAS_STARTING_OFFSET(n, incx);
	size_t iy = TBLAS_STARTING_OFFSET(n, incy);
	for(size_t i = 0; i < n; i++){
		const ScalarRealT xx = x[ix];
		const ScalarRealT yy = y[iy];
		x[ix] = h11 * xx + h12 * yy;
		y[iy] = h21 * xx + h22 * yy;
		ix += incx;
		iy += incy;
	}
}
