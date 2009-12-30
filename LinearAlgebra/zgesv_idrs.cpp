#include <cstddef>
#include <cmath>
#include <complex>
#include <iostream>

#define OUTPUT_MATLAB
//#define OUTPUT_MATHEMATICA

// Implementation of IDR(s) - Induced Dimension Reduction method

void printcomplex(const std::complex<double> &c){
	#ifdef OUTPUT_MATHEMATICA
	std::cout << c.real();
	if(c.imag() < 0){
		std::cout << "-I*" << -(c.imag());
	}else if(c.imag() > 0){
		std::cout << "+I*" << c.imag();
	}
	#endif
	#ifdef OUTPUT_MATLAB
	std::cout << c.real();
	if(c.imag() < 0){
		std::cout << " - " << -(c.imag());
		std::cout << "i";
	}else if(c.imag() > 0){
		std::cout << " + " << c.imag();
		std::cout << "i";
	}
	#endif
}
void print(size_t n, std::complex<double> *V){
	#ifdef OUTPUT_MATHEMATICA
	std::cout << "{";
	for(size_t i = 0; i < n; ++i){
		printcomplex(V[i]);
		if(i != n-1){
			std::cout << ", ";
		}
	}
	std::cout << "}";
	#endif
	#ifdef OUTPUT_MATLAB
	std::cout << "[";
	for(size_t i = 0; i < n; ++i){
		printcomplex(V[i]);
		if(i != n-1){
			std::cout << std::endl;
		}
	}
	std::cout << "]";
	#endif
}
void print(size_t m, size_t n, size_t lda, std::complex<double> *A){
	#ifdef OUTPUT_MATLAB
	std::cout << "[";
	for(size_t i = 0; i < m; ++i){
		std::cout << "";
		for(size_t j = 0; j < n; ++j){
			printcomplex(A[i+j*lda]);
			if(j != n-1){
				std::cout << ", ";
			}
		}
		std::cout << std::endl;
	}
	std::cout << "]";
	#endif
	#ifdef OUTPUT_MATHEMATICA
	std::cout << "{";
	for(size_t i = 0; i < m; ++i){
		std::cout << "{";
		for(size_t j = 0; j < n; ++j){
			printcomplex(A[i+j*lda]);
			if(j != n-1){
				std::cout << ", ";
			}
		}
		std::cout << "}";
		if(i != m-1){
			std::cout << "," << std::endl;
		}
	}
	std::cout << "}";
	#endif
}

double vecnorm(const size_t n, const std::complex<double>* v){
	double normv = 0;
	for(size_t i = 0; i < n; ++i){ normv += std::norm(v[i]); }
	return sqrt(normv);
}

void orth(const size_t n, const size_t s, std::complex<double> *Q){
	for(size_t j = 0; j < s; ++j){
		double inorm = double(1)/vecnorm(n, &Q[0+j*n]);
		for(size_t i = 0; i < n; ++i){
			Q[i+j*n] *= inorm;
		}
		for(size_t k = j+1; k < s; ++k){
			std::complex<double> dot = 0;
			for(size_t i = 0; i < n; ++i){
				dot += std::conj(Q[i+j*n]) * Q[i+k*n];
			}
			for(size_t i = 0; i < n; ++i){
				Q[i+k*n] -= dot*Q[i+j*n];
			}
		}
		//std::cout << "Q[" << j << "] = "; print(n, &Q[0+j*n]); std::cout << std::endl;
	}
}


size_t lu(const size_t m, const size_t n, const size_t lda, std::complex<double> *A, size_t *pivots){
	size_t info = 0;
	const size_t min_dim = (m<n) ? m:n;
	for(size_t j = 0; j < min_dim; ++j){
		size_t jp = j;
		{
			double largestval = std::abs(A[j+j*lda]);
			for(size_t i = j+1; i < m; ++i){
				double curval = std::abs(A[i+j*lda]);
				if(curval > largestval){
					jp = i;
					largestval = curval;
				}
			}
		}
		pivots[j] = jp;
		if(double(0) != A[jp+j*lda]){
			if(jp != j){
				for(size_t i = 0; i < n; ++i){
					std::swap(A[j+i*lda], A[jp+i*lda]);
				}
			}
			for(size_t i = j+1; i < m; ++i){
				A[i+j*lda] /= A[j+j*lda];
			}
		}else{
			info = j+1;
		}
		for(size_t q = j+1; q < m; ++q){
			for(size_t p = j+1; p < n; ++p){
				A[p+q*lda] -= A[p+j*lda]*A[j+q*lda];
			}
		}
	}
	return info;
}
void lu_solve(const size_t m, const size_t n, const size_t lda, std::complex<double> *A, size_t *pivots, std::complex<double> *b){
	// Apply pivots
	for(size_t i = 0; i < m; ++i){
		if(pivots[i] != i){
			std::swap(b[i], b[pivots[i]]);
		}
	}
	// Solve lower unit
	for(size_t k = 0; k < m; ++k){
		if(double(0) != b[k]){
			for(size_t i = k+1; i < m; ++i){
				b[i] -= b[k]*A[i+k*lda];
			}
		}
	}
	// Solve upper non unit
	for(int k = m-1; k >= 0; --k){
		if(double(0) != b[k]){
			b[k] /= A[k+k*lda];
			for(int i = 0; i < k; ++i){
				b[i] -= b[k]*A[i+k*lda];
			}
		}
	}
}
int zgesv_idrs(
	const size_t n,
	// A is a function which multiplies the matrix by the first argument
	// and returns the result in the second. The second argument must
	// be manually cleared. The third parameter is user data, passed in
	// through Adata.
	void (*A)(const std::complex<double>*, std::complex<double>*, void*),
	std::complex<double>* b,
	std::complex<double>* x,
	// Optional parameters
	void *Adata = NULL,
	size_t maxit = 0, // default is min(2*n,1000)
	const size_t s = 4,
	const double tol = 1e-8,
	bool x_initialized = false,
	// P is a precondition which simply solves P*x' = x,
	// where x i the first argument. The second parameter is user data,
	// which is passed in through Pdata.
	void (*P)(std::complex<double>*, void*) = NULL,
	void *Pdata = NULL,
	double angle = 0.7
){
	double normb = vecnorm(n, b);
	if(0 == normb){
		for(size_t i = 0; i < n; ++i){ x[i] = 0; }
		return 0;
	}
	const double tolb = tol*normb; // compute tolerance
	
	// Set initial x
	if(!x_initialized){
		for(size_t i = 0; i < n; ++i){ x[i] = 0; }
	}
	
	
	std::complex<double> *r = new std::complex<double>[n];
	A(x,r,Adata);
	for(size_t i = 0; i < n; ++i){ r[i] = b[i]-r[i]; }
	double normr = vecnorm(n, r);
	// Now, r = b-A*x
	
	std::complex<double> *Q = new std::complex<double>[n*s];
	{ // set up shadow space
		
		for(size_t j = 0; j < s; ++j){
			for(size_t i = 0; i < n; ++i){
				Q[i+j*n] = (double)rand()/(double)RAND_MAX - 0.5;
			}
		}
		// Orthogonalize Q
		orth(n, s, Q);
	}
	
	std::complex<double> *G = new std::complex<double>[n*s];
	std::complex<double> *U = new std::complex<double>[n*s];
	std::complex<double> *M = new std::complex<double>[s*s];
	std::complex<double> *Mcopy = new std::complex<double>[s*s];
	size_t *pivots = new size_t[s];
	for(size_t j = 0; j < s; ++j){
		for(size_t i = 0; i < n; ++i){
			G[i+j*n] = 0;
			U[i+j*n] = 0;
		}
		for(size_t i = 0; i < s; ++i){
			if(i == j){
				M[i+j*s] = 1;
			}else{
				M[i+j*s] = 0;
			}
		}
	}
	std::complex<double> *f = new std::complex<double>[s];
	std::complex<double> *c = new std::complex<double>[s];
	std::complex<double> *v = new std::complex<double>[n];
	std::complex<double> *t = new std::complex<double>[n];
	size_t iter = 0;
	std::complex<double> om = 1;
	
	if(0 == maxit){
		maxit = 2*n;
		if(1000 < maxit){ maxit = 1000; }
	}
	
	int ret = 0;
	while(normr > tolb && iter < maxit){
		std::cout << "iter = " << iter << std::endl;
		
		// generate RHS for small system
		for(size_t j = 0; j < s; ++j){
			std::complex<double> sum = 0;
			for(size_t i = 0; i < n; ++i){
				sum += r[i] * std::conj(Q[i+j*n]);
			}
			f[j] = sum;
		}
		
		for(size_t k = 0; k < s; ++k){
			// solve small systems of M(k:s,k:s)*c(k:s) = f(k:s)
			{
				// Copy over stuff for a destructive LU solve in Mcopy
				for(size_t j = k; j < s; ++j){
					for(size_t i = k; i < s; ++i){
						Mcopy[i+j*s] = M[i+j*s];
					}
					c[j] = f[j];
				}
				// Perform LU solve...
				lu(s-k, s-k, s, &Mcopy[k+k*s], pivots);
				lu_solve(s-k, s-k, s, &Mcopy[k+k*s], pivots, &c[k]);
			}
			// v = r - G(:,k:s)*c;
			for(size_t i = 0; i < n; ++i){
				std::complex<double> sum = 0;
				for(size_t j = k; j < s; ++j){
					sum += G[i+j*n]*c[j];
				}
				v[i] = r[i] - sum;
			}
			if(NULL != P){
				P(v, Pdata);
			}
			
			//U(:,k) = U(:,k:s)*c + om*v;
			for(size_t i = 0; i < n; ++i){
				std::complex<double> sum = 0;
				for(size_t j = k; j < s; ++j){
					sum += U[i+j*n]*c[j];
				}
				U[i+k*n] = sum + om*v[i];
			}
			//G(:,k) = A*U(:,k);
			A(&U[0+k*n], &G[0+k*n], Adata);
			
			// Bi-Orthogonalise the new basis vectors
			for(size_t j = 0; j < k; ++j){
				std::complex<double> alpha = 0;
				for(size_t i = 0; i < n; ++i){
					alpha += std::conj(Q[i+j*n])*G[i+k*n];
				}
				alpha /= M[j+j*s];
				for(size_t i = 0; i < n; ++i){
					G[i+k*n] -= alpha*G[i+j*n];
				}
				for(size_t i = 0; i < n; ++i){
					U[i+k*n] -= alpha*U[i+j*n];
				}
			}
			// New column of M = (Q'*G)'  (first k-1 entries are zero)
			for(size_t j = k; j < s; ++j){
				std::complex<double> sum = 0;
				for(size_t i = 0; i < n; ++i){
					sum += G[i+k*n]*std::conj(Q[i+j*n]);
				}
				M[j+k*s] = sum;
			}

			// Make r orthogonal to p_i, i = 1..k
			std::complex<double> beta = f[k]/M[k+k*s];
			for(size_t i = 0; i < n; ++i){
				r[i] -= beta*G[i+k*n];
			}
			for(size_t i = 0; i < n; ++i){
				x[i] += beta*U[i+k*n];
			}

			++iter;
			normr = vecnorm(n, r);

			if(normr < tolb || iter == maxit){ break; }
			
			// New f = Q'*r (first k  components are zero)
			for(size_t j = k+1; j < s; ++j){
				f[j] -= beta*M[j+k*s];
			}
		} // end k loop
		
		// If we break'd out of the inner loop, do so again
		if(normr < tolb){ break; }

		// Now we have sufficient vectors in G_j to compute residual in G_j+1
		// Note: r is already perpendicular to Q so v = r
		for(size_t i = 0; i < n; ++i){ v[i] = r[i]; }
		if(NULL != P){
			P(v, Pdata);
		}
		A(v, t, Adata);
		{ // compute new omega
			double norms = vecnorm(n, r), normt = vecnorm(n, t);
			std::complex<double> ts = 0;
			for(size_t i = 0; i < n; ++i){
				ts += std::conj(t[i])*r[i];
			}
			double rho = std::abs(ts/(normt*norms));
			om = ts/(normt*normt);
			if(rho < angle){
				om *= angle/rho;
			}
		}
		
		for(size_t i = 0; i < n; ++i){ r[i] -= om*t[i]; }
		for(size_t i = 0; i < n; ++i){ x[i] += om*v[i]; }
		normr = vecnorm(n, r);
		++iter;
	}
	
	delete [] r;
	delete [] G;
	delete [] U;
	delete [] M;
	delete [] Mcopy;
	delete [] f;
	delete [] c;
	delete [] v;
	delete [] t;
	return ret;
}

#include <map>
struct ccsmatrix{
	int cols;
	int *colptr;
	int *rowind;
	std::complex<double> *values;
	
	typedef std::pair<size_t,size_t> index_t;
	struct index_comp_t{ // sort by columns first, then rows
		bool operator()(const index_t &a, const index_t &b) const{
			if(a.second < b.second){ return true; }
			else if(a.second > b.second){ return false; }
			else{ return a.first < b.first; }
		}
	};
	typedef std::map<index_t,std::complex<double>,index_comp_t> entry_map_t;
	// Assumes that entries contains at least one element per column.
	// Assumes all indexes in entries is consistent with the r and c given.
	// Assumes that flags is set consistent with entries.
	ccsmatrix(size_t r, size_t c, const entry_map_t &entries):cols(c){
		size_t nnz = entries.size();
		values = new std::complex<double>[nnz];
		colptr = new int[cols+1];
		rowind = new int[nnz];
		
		size_t ip = 0;
		size_t prevcol = 0;
		colptr[0] = 0;
		for(ccsmatrix::entry_map_t::const_iterator i = entries.begin(); i != entries.end(); ++i){
			size_t col = i->first.second;
			rowind[ip] = i->first.first;
			values[ip] = i->second;
			++ip;
			if(prevcol != col){
				prevcol = col;
				colptr[col] = ip-1;
			}
		}
		colptr[cols] = nnz; // do this at the end in case entries was bad, at least this is correct
	}
	~ccsmatrix(){
		if(NULL != values){ delete [] values; }
		if(NULL != rowind){ delete [] rowind; }
		if(NULL != colptr){ delete [] colptr; }
	}
};

inline double frand(){
	return (double)rand()/(double)RAND_MAX - 0.5;
}
inline std::complex<double> crand(){
	return std::complex<double>(frand(),frand());
}
void print(const ccsmatrix *A){
	#ifdef OUTPUT_MATHEMATICA
	std::cout << "{";
	for(int j = 0; j < (int)A->cols; j++){
		for(int ip = A->colptr[j]; ip < A->colptr[j+1]; ip++){
			int i = A->rowind[ip];
			std::cout << "{" << i+1 << ", " << j+1 << "} -> ";
			printcomplex(A->values[ip]);
			std::cout << ", ";
		}
	}
	std::cout << "{_, _} -> 0}";
	#endif
	#ifdef OUTPUT_MATLAB
	std::cout << "spconvert([";
	for(int j = 0; j < (int)A->cols; j++){
		for(int ip = A->colptr[j]; ip < A->colptr[j+1]; ip++){
			int i = A->rowind[ip];
			std::cout << i+1 << "\t" << j+1 << "\t";
			printcomplex(A->values[ip]);
			std::cout << std::endl;
		}
	}
	std::cout << "])";
	#endif
}

void mult(const ccsmatrix *A, const std::complex<double> *x, std::complex<double> *y){
	for(int i = 0; i < A->cols; ++i){
		y[i] = 0;
	}
	for(int j = 0; j < (int)A->cols; j++){
		for(int ip = A->colptr[j]; ip < A->colptr[j+1]; ip++){
			int i = A->rowind[ip];
			y[i] += x[j]*A->values[ip];
		}
	}
}

void Aop(const std::complex<double> *x, std::complex<double> *y, void *data){ 
	const ccsmatrix *A = reinterpret_cast<ccsmatrix*>(data);
	mult(A, x, y);
}


int main(){
	srand(time(0));
	size_t n = 100;
	std::complex<double> *b = new std::complex<double>[n];
	std::complex<double> *x = new std::complex<double>[n];
	std::complex<double> *r = new std::complex<double>[n];
	ccsmatrix::entry_map_t entries;
	for(size_t i = 0; i < n; ++i){
		entries[ccsmatrix::index_t(i,i)] = 1;
		b[i] = crand();
	}
	
	for(size_t i = 0; i < n; ++i){
		int row = rand()%n;
		int col = rand()%n;
		entries[ccsmatrix::index_t(row,col)] += crand();
	}
	ccsmatrix A(n, n, entries);

	//print(&A); std::cout << std::endl << std::endl;
	//print(n, b); std::cout << std::endl;
	
	zgesv_idrs(n, &Aop, b, x, (void*)&A, 100, 2);
	
	//print(n, x); std::cout << std::endl;
	
	mult(&A, x, r);
	for(size_t i = 0; i < n; ++i){
		r[i] -= b[i];
	}
	double normr = vecnorm(n, r);
	double normb = vecnorm(n, b);
	std::cout << "normb = " << normb << std::endl;
	std::cout << "err = " << normr/normb << std::endl;
	
	delete [] b;
	delete [] x;
	delete [] r;
	return 0;
}
