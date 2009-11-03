template <typename TypeAlpha, typename TypeA, typename TypeB,
          typename TypeBeta, typename TypeC>
void TBLAS_NAME(gemm,GEMM)(
  const char *transa, const char *transb,
  const size_t &m, const size_t &n, const size_t &k,
  const TypeAlpha &alpha,
  const TypeA *A, const size_t &lda,
  const TypeB *B, const size_t &ldb,
  const TypeBeta &beta,
  const TypeC *C, const size_t &ldc){

	if(0 == m || 0 == n || 0 == k){ return; }
	if(0 == alpha && 1 == beta){ return; }

	unsigned int bits = 0x3;
	if(*transa == TBLAS::Op::None){
		bits &= 0xE;
	}else if(*transa == TBLAS::Op::ConjugateTranspose){ bits |= 0x4; }
	if(*transb == TBLAS::Op::None){
		bits &= 0xD;
	}else if(*transb == TBLAS::Op::ConjugateTranspose){ bits |= 0x8; }
	
	switch(bits){
	case 0x0: // N, N
		for(size_t i = 0; i < m; ++i){
			for(size_t j = 0; j < n; ++j){
				TC::value_type sum(0);
				for(size_t q = 0; q < k; ++q){
					sum += A(i,q)*B(q,j);
				}
				C(i,j) = scale_AB * sum + scale_C * C(i,j);
			}
		}
		break;
	case 0x1: // T, N
		for(size_t i = 0; i < m; ++i){
			for(size_t j = 0; j < n; ++j){
				TC::value_type sum(0);
				for(size_t q = 0; q < k; ++q){
					sum += A(q,i)*B(q,j);
				}
				C(i,j) = scale_AB * sum + scale_C * C(i,j);
			}
		}
		break;
	case 0x2: // N, T
		for(size_t i = 0; i < m; ++i){
			for(size_t j = 0; j < n; ++j){
				TC::value_type sum(0);
				for(size_t q = 0; q < k; ++q){
					sum += A(i,q)*B(j,q);
				}
				C(i,j) = scale_AB * sum + scale_C * C(i,j);
			}
		}
		break;
	case 0x3: // T, T
		for(size_t i = 0; i < m; ++i){
			for(size_t j = 0; j < n; ++j){
				TC::value_type sum(0);
				for(size_t q = 0; q < k; ++q){
					sum += A(q,i)*B(j,q);
				}
				C(i,j) = scale_AB * sum + scale_C * C(i,j);
			}
		}
		break;
	case 0x5: // C, N
		for(size_t i = 0; i < m; ++i){
			for(size_t j = 0; j < n; ++j){
				TC::value_type sum(0);
				for(size_t q = 0; q < k; ++q){
					sum += TBLAS_CONJ(A(q,i))*B(q,j);
				}
				C(i,j) = scale_AB * sum + scale_C * C(i,j);
			}
		}
		break;
	case 0xA: // N, C
		for(size_t i = 0; i < m; ++i){
			for(size_t j = 0; j < n; ++j){
				TC::value_type sum(0);
				for(size_t q = 0; q < k; ++q){
					sum += A(i,q)*TBLAS_CONJ(B(j,q));
				}
				C(i,j) = scale_AB * sum + scale_C * C(i,j);
			}
		}
		break;
	case 0xF: // C, C
		for(size_t i = 0; i < m; ++i){
			for(size_t j = 0; j < n; ++j){
				TC::value_type sum(0);
				for(size_t q = 0; q < k; ++q){
					sum += A(q,i)*B(j,q);
				}
				C(i,j) = scale_AB * TBLAS_CONJ(sum) + scale_C * C(i,j);
			}
		}
		break;
	default:
		TBLAS_ERROR("Unrecognized operation");
		break;
	}
}
