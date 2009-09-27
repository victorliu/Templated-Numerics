template <typename ScalarT>
void TBLAS_NAME(copy,COPY)(
  const size_t &n,
  const ScalarT *src, const size_t &incsrc,
  ScalarT *dst, const size_t &incdst){
	size_t ix = TBLAS_STARTING_OFFSET(n, incsrc);
	size_t iy = TBLAS_STARTING_OFFSET(n, incdst);
	for(size_t i = 0; i < n; ++i){
		dst[iy] = src[ix];
		ix += incsrc;
		iy += incdst;
	}
}
