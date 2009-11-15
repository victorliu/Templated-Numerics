#ifndef _TBLAS_DEFS_H_
#define _TBLAS_DEFS_H_

#define TBLAS_UINT unsigned int
#define TBLAS_INT int

#define TBLAS_EXTENSIONS // want to use the extra functions
#define TBLAS_OVERFLOW_PROTECTION // want to go prevent against over/underflow

#define TBLAS_NAME(lower,upper) lower
#define TBLAS_ASSERT(COND, DESC) do{ if(!(COND)){ TBLAS_ERROR(DESC); } }while(0)
#define TBLAS_ERROR(DESC) fprintf(stderr, "Error in %s: %s", __FUNCTION__, DESC)

#define TBLAS_STARTING_OFFSET(N, INC) ((INC)>0 ? 0 : ((N) - 1) * (-(INC)))
#define GB(KU,KL,lda,i,j) ((KU+1+(i-j))*lda + j)
#define TRCOUNT(N,i) ((((i)+1)*(2*(N)-(i)))/2)
/* #define TBUP(N,i,j) */
/* #define TBLO(N,i,j) */
#define TPUP(N,i,j) (TRCOUNT(N,(i)-1)+(j)-(i))
#define TPLO(N,i,j) (((i)*((i)+1))/2 + (j))


#define TBLAS_ABS(arg) std::abs(arg)
#define TBLAS_SQRT(arg) std::sqrt(arg)
#define TBLAS_CONJ(arg) std::conj(arg)

#define TBLAS_TRAITS(type) TypeTraits<type>

#endif
