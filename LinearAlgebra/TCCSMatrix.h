#ifndef _TCCSMATRIX_H_
#define _TCCSMATRIX_H_

#define USING_TCCS_MATRIX

template <typename NumericType, class TAllocator = std::allocator<NumericType> >
class TCCSMatrix : public ReadableMatrix<NumericType>{
public:
	// Compatible with TAUCS
	int cols, rows; // rows is unused if symmetric
	int flags;
	int *colptr;
	int *rowind;
	NumericType *values;
	TAllocator allocator;

	typedef NumericType value_type;
	
	typedef ReadableMatrix<value_type> readable_matrix;
	typedef WritableMatrix<value_type> non_writable_matrix;
	
	// flag values
	static const int LOWER = 1;
	static const int UPPER = 2;
	static const int TRIANGULAR = 4;
	static const int SYMMETRIC = 8;
	static const int HERMITIAN = 16;
	static const int PATTERN = 32;
	
	TCCSMatrix():cols(0),rows(0),flags(0),colptr(NULL),rowind(NULL),values(NULL){}
	TCCSMatrix(size_t r, size_t c, size_t nnz):cols(c),rows(r),flags(0),colptr(NULL),rowind(NULL),values(NULL){
		values = allocator.allocate(nnz);
		colptr = new int[cols+1];
		rowind = new int[nnz];
	}
	TCCSMatrix(const TCCSMatrix &M):cols(M.Cols()),rows(M.Rows()),flags(M.Flags()),colptr(NULL),rowind(NULL),values(NULL){
		values = allocator.allocate(M.nnz());
		colptr = new int[cols+1];
		rowind = new int[M.nnz()];
		
		std::uninitialized_copy(M.colptr, M.colptr+M.Cols()+1, colptr);
		std::uninitialized_copy(M.rowind, M.rowind+M.nnz(), rowind);
		std::uninitialized_copy(M.values, M.values+M.nnz(), values);
	}
	
	typedef std::pair<size_t,size_t> index_t;
	struct index_comp_t{ // sort by columns first, then rows
		void operator()(const index_t &a, const index_t &b) const{
			if(a.second < b.second){ return true; }
			else if(a.second > b.second){ return false; }
			else{ return a.first < b.first; }
		}
	};
	typedef std::map<index_t,value_type,index_comp_t> entry_map_t;
	// Assumes that entries contains at least one element per column.
	// Assumes all indexes in entries is consistent with the r and c given.
	// Assumes that flags is set consistent with entries.
	TCCSMatrix(size_t r, size_t c, const entry_map_t &entries, int flags):cols(c),rows(r),flags(flags){
		size_t nnz = entries.size();
		values = allocator.allocate(nnz);
		colptr = new int[cols+1];
		rowind = new int[nnz];
		
		size_t ip = 0;
		int prevcol = 0;
		colptr[0] = 0;
		for(entry_map_t::const_iterator i = entries.begin(); i != entries.end(); ++i){
			int col = i->first.second;
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
	TCCSMatrix& operator=(const TCCSMatrix &M){
		if(NULL != rowind){ delete [] rowind; }
		if(NULL != values){ allocator.deallocate(values, nnz()); }
		if(NULL != colptr){ delete [] colptr; }

		cols = M.Cols();
		rows = M.Rows();
		colptr = new int[cols+1];
		std::uninitialized_copy(M.colptr, M.colptr+M.cols+1, colptr);

		values = allocator.allocate(nnz());
		rowind = new int[nnz()];
		
		std::uninitialized_copy(M.rowind, M.rowind+M.nnz(), rowind);
		std::uninitialized_copy(M.values, M.values+M.nnz(), values);

		return *this;
	}
	virtual ~TCCSMatrix(){
		if(NULL != values){ allocator.deallocate(values, nnz()); }
		if(NULL != rowind){ delete [] rowind; }
		if(NULL != colptr){ delete [] colptr; }
	}
	
	size_t Rows() const{ return rows; }
	size_t Cols() const{ return cols; }
	size_t nnz() const{ return colptr[cols]; }
	
	value_type operator()(size_t row, size_t col) const{
#ifdef USE_MATRIX_OPS_ASSERTS
		assert(row < Rows());
		assert(col < Cols());
#endif
		int col_start = colptr[col];
		for(int i = col_start; i < colptr[col+1]; ++i){
			if((int)row == rowind[i]){
				return values[i];
			}else if((int)row > rowind[i]){
				break;
			}
		}
		return value_type(0);
	}
};

#endif // _TCCSMATRIX_H_
