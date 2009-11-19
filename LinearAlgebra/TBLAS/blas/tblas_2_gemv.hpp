template <typename TypeAlpha, typename TypeA, typename TypeX,
          typename TypeBeta, typename TypeY>
void TBLAS_NAME(gemv,GEMV)(
  const char *trans,
  const size_t &m, const size_t &n,
  const TypeAlpha& alpha,
  const TypeA *A, const size_t &lda,
  const TypeX *x, const size_t &incx,
  const TypeBeta& beta,
  const TypeY *y, const size_t &incy){
	
	if(0 == m || 0 == n){ return; }
	if(0 == alpha && 0 == beta){ return; }
	
	size_t lenx, leny;
	if(trans[0] == TBLAS::Op::None){
		lenx = n;
		leny = m;
	}else{
		lenx = m;
		leny = n;
	}

	// y <- beta*y
	if(beta == 0){
		size_t iy = TBLAS_STARTING_OFFSET(leny, incy);
		for(size_t i = 0; i < leny; i++){
			y[iy] = 0;
			iy += incy;
		}
	}else if(beta != 1){
		size_t iy = TBLAS_STARTING_OFFSET(leny, incy);
		for(size_t i = 0; i < leny; i++){
			y[iy] *= beta;
			iy += incy;
		}
	}

	if(alpha == 0.0){ return; }

	if(trans[0] == TBLAS::Op::None){
		// y <- alpha*A*x + y
		size_t ix = TBLAS_STARTING_OFFSET(lenx, incx);
		for (size_t j = 0; j < lenx; j++) {
			const TypeX temp = alpha * x[ix];
			if(temp != 0){
				size_t iy = TBLAS_STARTING_OFFSET(leny, incy);
				for(size_t i = 0; i < leny; i++) {
					y[iy] += temp * A[lda * j + i];
					iy += incy;
				}
			}
			ix += incx;
		}
	}else if(trans[0] == TBLAS::Op::Transpose){
		// y <- alpha*A'*x + y
		size_t iy = TBLAS_STARTING_OFFSET(leny, incy);
		for(size_t i = 0; i < leny; i++){
			TypeY temp = 0;
			size_t ix = TBLAS_STARTING_OFFSET(lenx, incx);
			for(size_t j = 0; j < lenx; j++){
				temp += x[ix] * A[lda*i + j];
				ix += incx;
			}
			y[iy] += alpha * temp;
			iy += incy;
		}
	}else if(trans[0] == TBLAS::Op::ConjugateTranspose){
		// y <- alpha*A^H*x + y
		size_t iy = TBLAS_STARTING_OFFSET(leny, incy);
		for(size_t i = 0; i < leny; i++){
			TypeX temp = 0;
			size_t ix = TBLAS_STARTING_OFFSET(lenx, incx);
			for(size_t j = 0; j < lenx; j++){
				temp += x[ix] * TBLAS_CONJ(A[lda*i + j]);
				ix += incx;
			}
			y[iy] += alpha * temp;
			iy += incy;
		}
	}else{
		TBLAS_ERROR("Unrecognized operation");
	}
}
