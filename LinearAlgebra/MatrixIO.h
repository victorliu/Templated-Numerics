#ifndef _MATRIX_IO_H_
#define _MATRIX_IO_H_

#include <iostream>
#include "MatrixInterfaces.h"
#include "MatrixViews.h"

// Dense Mathematica output
template <class T>
	typename IsReadableMatrix<typename T::readable_matrix,
std::ostream&
	>::type
operator<<(std::ostream& os, const T &view){
	os << "{" << std::endl;
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
	os << "}";
	return os;
}

// Dense Mathematica output
template <class T>
	typename IsReadableVector<typename T::readable_vector,
std::ostream&
	>::type
operator<<(std::ostream& os, const T &view){
	os << "{";
	const size_t N = view.size();
	for(size_t i = 0; ; ++i){
		os << view[i];
		if(N-1 == i){
			os << "}";
			break;
		}else{
			os << ", ";
		}
	}
	return os;
}

#endif // _MATRIX_IO_H_
