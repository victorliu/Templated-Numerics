#ifndef _MATRIX_IO_H_
#define _MATRIX_IO_H_

namespace MatrixIO{

// Dense Mathematica output
std::ostream& operator<<(std::ostream& os, const MatrixViewBase &view){
	os << "{" std::endl;
	const size_t M = view.Rows(), N = view.Cols();
	for(size_t i = 0; i < M; ++i){
		os << "{";
		for(size_t j = 0; ; ++j){
			os << view(i,j);
			if(N-1 == j){
				if(i < M-1){
					os << "}," << std::endl;
				}else{
					os << "}" << std::endl;
				}
				break;
			}else{
				os << ", ";
			}
		}
	}
	os << "}" std::endl;
	return *os;
}

// Dense Mathematica output
std::ostream& operator<<(std::ostream& os, const VectorViewBase &view){
	os << "{";
	const size_t N = view.Cols();
	for(size_t i = 0; ; ++i){
		os << view[i];
		if(N-1 == i){
			os << "}" << std::endl;
			break;
		}else{
			os << ", ";
		}
	}
	return *os;
}

}; // namespace MatrixIO

#endif // _MATRIX_IO_H_
