#ifndef _TONES_VECTOR_H_
#define _TONES_VECTOR_H_

#define USING_TONES_VECTOR

#include "MatrixInterfaces.h"

template <typename NumericType>
class TOnesVector : public ReadableVector<NumericType>{
	size_t rows;
	TOnesVector& operator=(const TOnesVector &V){ return *this; }
public:
	typedef NumericType value_type;
	typedef ReadableVector<value_type> non_view_type;
	typedef ReadableVector<value_type> readable_vector;
	typedef WritableVector<value_type> non_writable_vector;
	
	TOnesVector():rows(0){}
	TOnesVector(size_t r):rows(r){}
	TOnesVector(const TOnesVector &V):rows(V.rows){}
	virtual ~TOnesVector(){}
	
	void Resize(size_t nRows){ rows = nRows; }
	
	size_t size() const{ return rows; }
	value_type operator[](size_t row) const{ return value_type(1); }
};

#endif // _TONES_VECTOR_H_
