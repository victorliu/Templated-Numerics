#ifndef _ALLOCATOR_ADAPTER_H_
#define _ALLOCATOR_ADAPTER_H_

namespace std{

template <class T, class MallocFunc, class FreeFunc>
class allocator_adapter{
	static size_t allocated;
	static size_t max_allocated;
	MallocFunc m;
	FreeFunc f;
public:
	// type definitions
	typedef T        value_type;
	typedef T*       pointer;
	typedef const T* const_pointer;
	typedef T&       reference;
	typedef const T& const_reference;
	typedef std::size_t    size_type;
	typedef std::ptrdiff_t difference_type;
	//static size_t allocated;

	// rebind allocator to type U
	template <class U>
	struct rebind{
		typedef allocator_adapter<U, MallocFunc, FreeFunc> other;
	};

	// return address of values
	pointer address(reference value) const{ return &value; }
	const_pointer address(const_reference value) const{ return &value; }

	allocator_adapter() throw(){}
	allocator_adapter(const allocator_adapter&) throw(){}
	template <class U>
	allocator_adapter(const allocator_adapter<U,MallocFunc,FreeFunc>&) throw(){}
	~allocator_adapter() throw(){}

	size_type max_size() const throw(){ return std::numeric_limits<std::size_t>::max() / sizeof(T); }

	// allocate but don't initialize num elements of type T
	pointer allocate(size_type num, const void* = 0){
		pointer ret = (pointer)m(num*sizeof(T)); // (pointer)(::operator new(num*sizeof(T)));
		allocated += num * sizeof(T);
		if(allocated > max_allocated){ max_allocated = allocated; }
		return ret;
	}

	// initialize elements of allocated storage p with value value
	void construct(pointer p, const T& value){
		new((void*)p)T(value);
	}

	void destroy(pointer p){ p->~T(); }

	// deallocate storage p of deleted elements
	void deallocate(pointer p, size_type num){
		f((void*)p); // ::operator delete((void*)p);
		allocated -= num * sizeof(T);
	}
};

template <class T, class MallocFunc, class FreeFunc>
size_t allocator_adapter<T,MallocFunc,FreeFunc>::allocated = 0;

template <class T, class MallocFunc, class FreeFunc>
size_t allocator_adapter<T,MallocFunc,FreeFunc>::max_allocated = 0;

template <class MallocFunc, class FreeFunc>
class allocator_adapter<void,MallocFunc,FreeFunc>{
public:
	typedef size_t      size_type;
	typedef ptrdiff_t   difference_type;
	typedef void*       pointer;
	typedef const void* const_pointer;
	typedef void        value_type;

	template<typename _Tp1>
	struct rebind{ typedef allocator_adapter<_Tp1,MallocFunc,FreeFunc> other; };
};

template <class T1, class T2, class MallocFunc, class FreeFunc>
bool operator== (const allocator_adapter<T1,MallocFunc,FreeFunc>&, const allocator_adapter<T2,MallocFunc,FreeFunc>&) throw(){ return true; }
template <class T1, class T2, class MallocFunc, class FreeFunc>
bool operator!= (const allocator_adapter<T1,MallocFunc,FreeFunc>&, const allocator_adapter<T2,MallocFunc,FreeFunc>&) throw(){ return false; }

}; // namespace std

#endif // _ALLOCATOR_ADAPTER_H_
