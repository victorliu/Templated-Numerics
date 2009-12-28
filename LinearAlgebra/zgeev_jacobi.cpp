#include <cstddef>
#include <cmath>
#include <complex>
#include <iostream>

// A non-parallel implementation of:
//   A parallel algorithm for the eigenvalues
//   and eigenvectors of a general complex matrix
//   by Guatam M. Shroff
//   Numerische Mathematik 58, pages 779-805 (1991)
// The matrix is only reduced to upper triangular form,
// but in practice it is always diagonal.

// Returns the unitary transformation parameters
// Returns the reduction in lower triangular norm
static inline double
pnrj_unitary(size_t n, size_t lda, std::complex<double> *A, size_t p, size_t q, std::complex<double> cs[3]){
	std::complex<double> Apq = A[p+q*lda];
	double aApq2 = std::norm(Apq);
	if(0 != aApq2){
		std::complex<double> Aqp = A[q+p*lda];
		std::complex<double> dpq = 0.5*(A[p+p*lda] - A[q+q*lda]);
		std::complex<double> root = sqrt(dpq*dpq + Apq*Aqp);
		
		std::complex<double> dmax_plus  = dpq + root;
		std::complex<double> dmax_minus = dpq - root;
		std::complex<double> dmax;
		if(std::norm(dmax_plus) > std::norm(dmax_minus)){
			dmax = dmax_plus;
		}else{
			dmax = dmax_minus;
		}
		// At this point we should choose
		// tan x = exp(i theta) Apq/dmax
		// such that theta makes tan x real
		// The resulting parameters are then
		// [ cs[0]  cs[2] ]
		// [ cs[1]  cs[3] ]
		// where
		// cs[0] = cs[3] = cos x, cs[1] = exp(i theta) sin x
		// cs[2] = -exp(-i theta) sin x
		//
		// Thus, tan x = std::abs(Apq/dmax)
		dmax = Apq/dmax; // we use dmax as Aqp/dmax
		double tanx = std::abs(dmax);
		std::complex<double> exp_i_theta = tanx/dmax; // this is actually exp(-i theta)
		if(tanx > 1){ tanx = 1; } // limit the maximum rotation angle
		// Now cos^2 = 1/(1+tan^2), sin^2 = tan^2 cos^2
		cs[0] = double(1)/sqrt(tanx*tanx+1);
		//cs[3] = cs[0];
		cs[1] = tanx*cs[0]; // cs[1] = sin x
		cs[2] = cs[1];
		cs[1] *= exp_i_theta;
		cs[2] /= -exp_i_theta;
	}else{
		cs[0] = 1; // cs[3] = cs[0];
		cs[1] = cs[2] = 0;
	}

	return aApq2;
}

static inline double
pnrj_shear(size_t n, size_t lda, std::complex<double> *A, size_t p, size_t q, std::complex<double> cs[3]){
	double Gpq = 0;
	std::complex<double> cpq(0);
	for(size_t j = 0; j < n; ++j){
		cpq += (A[p+j*lda]*std::conj(A[q+j*lda]) - std::conj(A[j+p*lda])*A[j+q*lda]);
		
		if(j == p || j == q){ continue; }
		double Gterm = 0;
		Gterm += std::norm(A[p+j*lda]);
		Gterm += std::norm(A[q+j*lda]);
		Gterm += std::norm(A[j+p*lda]);
		Gterm += std::norm(A[j+q*lda]);
		Gpq += Gterm;
	}
	std::complex<double> dpq = A[q+q*lda] - A[p+p*lda];
	// xi_pq = exp(i alpha) Aqp + exp(-i alpha) Apq
	// alpha = arg(cpq) - pi/2
	// Thus, xi_pq = -i exp(i arg(cpq)) Aqp + i exp(-i arg(cpq)) Apq
	// But exp(i arg(cpq)) is simply cpq/|cpq|
	double acpq = std::abs(cpq);
	if(0 == acpq){
		cs[0] = 1;
		cs[1] = cs[2] = 0;
		return 0;
	}
	std::complex<double> eialpha = std::complex<double>(0,-1)*(cpq/acpq);
	std::complex<double> xipq = eialpha*A[q+p*lda] + A[p+q*lda]/eialpha;
	// Now, we will generate the transformation
	// [ cs[0]  cs[2] ]
	// [ cs[1]  cs[3] ]
	// where
	// cs[0] = cs[3] = cosh y,
	// cs[1] =  i exp(-i alpha) sinh y
	// cs[2] = -i exp( i alpha) sinh y
	// and
	// tanh y = -|cpq| / (2*(|dpq|^2 + |xipq|^2) + Gpq)
	double tanhy = -acpq / (2*(std::norm(dpq) + std::norm(xipq)) + Gpq);
	// cosh^2 - sinh^2 = 1, tanh = sinh/cosh
	double coshy = double(1)/sqrt(double(1) - tanhy*tanhy);
	cs[0] = coshy; // cs[3] = cs[0]
	double sinhy = coshy*tanhy;
	cs[1] = std::complex<double>(0, sinhy)/eialpha;
	cs[2] = std::complex<double>(0,-sinhy)*eialpha;

	return 0;
}

static inline double
pnrj_diagonal(size_t n, size_t lda, std::complex<double> *A, size_t j, double *t){
	double g = 0;
	double h = 0;
	for(size_t i = 0; i < n; ++i){
		if(i == j){ continue; }
		g += std::norm(A[i+j*lda]);
		h += std::norm(A[j+i*lda]);
	}
	g = sqrt(g); h = sqrt(h);
	*t = sqrt(h/g);
	
	static const double tlimit = 1e8;
	
	if(*t > tlimit){
		*t = tlimit;
		h = tlimit*tlimit*g;
	}else if(*t < double(1)/tlimit){
		*t = double(1)/tlimit;
		g = h/(tlimit*tlimit);
	}
	
	// Compute the estimate of norm reduction
	g -= h;
	return g*g;
}

inline static void
pnrj_apply_rotation_L(size_t n, size_t lda, std::complex<double>* A, size_t p, size_t q, std::complex<double> cs[3]){
	// Apply rotation to matrix A,  A' = J^{-1} A
	for (size_t j = 0; j < n; j++){
		std::complex<double> Apj = A[p+j*lda];
		std::complex<double> Aqj = A[q+j*lda];
		A[p+j*lda] =  Apj * cs[0] - Aqj * cs[2];
		A[q+j*lda] = -Apj * cs[1] + Aqj * cs[0]; // cs[3]
	}
}

inline static void
pnrj_apply_rotation_R(size_t n, size_t lda, std::complex<double>* A, size_t p, size_t q, std::complex<double> cs[3]){
	// Apply rotation to matrix A,  A' = A J
	for(size_t i = 0; i < n; i++){
		std::complex<double> Aip = A[i+p*lda];
		std::complex<double> Aiq = A[i+q*lda];
		A[i+p*lda] = Aip * cs[0] + Aiq * cs[1];
		A[i+q*lda] = Aip * cs[2] + Aiq * cs[0]; // cs[3]
	}
}

inline static void
pnrj_diagonal_L(size_t n, size_t lda, std::complex<double>* A, size_t j, double t){
	t = (double)1/t;
	// Apply diagonal to matrix A,  A' = A J
	for(size_t i = 0; i < n; i++){
		A[j+i*lda] *= t;
	}
}

inline static void
pnrj_diagonal_R(size_t n, size_t lda, std::complex<double>* A, size_t j, double t){
	// Apply diagonal to matrix A,  A' = A J
	for(size_t i = 0; i < n; i++){
		A[i+j*lda] *= t;
	}
}

int zgeev_jacobi(size_t n, size_t lda, std::complex<double>* A, std::complex<double>* eval, size_t ldvec, std::complex<double>* evec, unsigned int max_iter, unsigned int *nrot){
	for(size_t p = 0; p < n; ++p){
		for(size_t q = 0; q < n; ++q){
			if(p == q){
				evec[q+p*ldvec] = 1;
			}else{
				evec[q+p*ldvec] = 0;
			}
		}
		eval[p] = 0;
	}

	double normL = std::numeric_limits<double>::max();
	double normL_prev;
	const double normL_threshold = n*n/2 * std::numeric_limits<double>::epsilon();

	size_t iter = 0;
	do{
		normL_prev = normL;
		// Compute the Frobenius norm of the lower triangle
		normL = 0;
		for(size_t q = 0; q < n-1; ++q){ // column
			for(size_t p = q+1; p < n; ++p){
				double absLpq2 = std::norm(A[p+q*lda]);
				normL += absLpq2;
			}
		}
		normL = sqrt(normL);
		std::cout << "norm = " << normL << std::endl;
		
		if(normL < normL_threshold){
			break;
		}
		
		// Perform a sweep
		//   A sweep is a set of rotations followed by a set of diagonal
		//   transformations. A rotation is a shear followed by a unitary
		//   transformation. Rotations are performed on all subdiagonal
		//   elements, while diagonal transformations are applied to each
		//   element of the diagonal.
		
		// Apply all rotations
		for(size_t q = 0; q < n-1; ++q){
			for(size_t p = q+1; p < n; ++p){
				std::complex<double> cs[3];
				
				pnrj_shear(n, lda, A, p, q, cs);
				pnrj_apply_rotation_L(n, lda, A, p, q, cs);
				pnrj_apply_rotation_R(n, lda, A, p, q, cs);
				pnrj_apply_rotation_R(n, ldvec, evec, p, q, cs);
				
				pnrj_unitary(n, lda, A, p, q, cs);
				pnrj_apply_rotation_L(n, lda, A, p, q, cs);
				pnrj_apply_rotation_R(n, lda, A, p, q, cs);
				pnrj_apply_rotation_R(n, ldvec, evec, p, q, cs);
			}
		}
		// Apply all diagonal transformations
		for(size_t j = 0; j < n; ++j){
			double t;
			pnrj_diagonal(n, lda, A, j, &t);
			pnrj_diagonal_L(n, lda, A, j, t);
			pnrj_diagonal_R(n, lda, A, j, t);
			pnrj_diagonal_R(n, ldvec, evec, j, t);
		}
		// max_iter ~<~ 2 + 2.8 log2(n)
	}while(iter++ <= max_iter);
	*nrot = iter;

	for(size_t p = 0; p < n; p++){
		eval[p] = A[p+p*lda];
	}

	if(iter == max_iter){ return 1; }

	return 0;
}

void print(size_t rows, size_t lda, std::complex<double> *A){
	std::cout << "{";
	for(size_t i = 0; i < rows; ++i){
		std::cout << "{";
		for(size_t j = 0; j < rows; ++j){
			std::cout << A[i+j*lda];
			if(j != rows-1){
				std::cout << ", ";
			}
		}
		std::cout << "}";
		if(i != rows-1){
			std::cout << ",";
		}
		std::cout << std::endl;
	}
	std::cout << "}";
}
void print(size_t n, std::complex<double> *V){
	std::cout << "{";
	for(size_t i = 0; i < n; ++i){
		std::cout << V[i];
		if(i != n-1){
			std::cout << ", ";
		}
	}
	std::cout << "}";
}



void makesym(size_t rows, size_t lda, std::complex<double> *A){
	for(size_t i = 0; i < rows; ++i){
		for(size_t j = i+1; j < rows; ++j){
			std::complex<double> v = 0.5*(A[i+j*lda]+A[j+i*lda]);
			A[i+j*lda] = A[j+i*lda] = v;
		}
	}
}

void makeid(size_t rows, size_t lda, std::complex<double> *A){
	for(size_t i = 0; i < rows; ++i){
		for(size_t j = 0; j < rows; ++j){
			if(i == j){
				A[j+i*lda] = 1;
			}else{
				A[j+i*lda] = 0;
			}
		}
	}
}

double frand(){
	return (double)rand()/(double)RAND_MAX;
}


int main(){
	srand(time(0));
	size_t n = 4;
	std::complex<double> *A = new std::complex<double>[n*n];
	
	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			A[i+j*n] = std::complex<double>(frand(), frand());
		}
	}
	/*
	makeid(n,n,A);
	A[1+0*n] = frand();
	A[2+1*n] = frand();
	A[2+0*n] = frand();
	makesym(n,n,A);
	*/
	
	//print(n,n,A); std::cout << std::endl;
	
	std::complex<double> *evals = new std::complex<double>[n];
	std::complex<double> *evecs = new std::complex<double>[n*n];
	size_t nrot;
	size_t max_iter = (int)(2+2.8*log((double)n)/log(2.0));
	if(max_iter < 8){ max_iter = 8; }
	std::cout << "theoretical num sweeps = " << max_iter << std::endl;
	zgeev_jacobi(n,n,A,evals,n,evecs, max_iter, &nrot);
	
	print(n,n,A); std::cout << std::endl;
	//print(n,evals); std::cout << std::endl;
	//print(n,n,evecs); std::cout << std::endl;
	
	delete [] A;
	delete [] evecs;
	delete [] evals;
	return 0;
}



