#ifndef _MEMORY_ALIGNMENT_H_
#define _MEMORY_ALIGNMENT_H_

template <unsigned int nByteAlign> // nByteAlign must be a power of 2
struct AlignedAllocate{
	void* operator()(size_t size){
		static const size_t align_size = nByteAlign;
		static const size_t align_mask = align_size - 1;

		char *ptr = (char *)malloc(size + align_size + sizeof(size_t));
		if(NULL == ptr){ return(NULL); }

		char *ptr2 = ptr + sizeof(size_t);
		char *aligned_ptr = ptr2 + (align_size - ((size_t)ptr2 & align_mask));

		ptr2 = aligned_ptr - sizeof(size_t);
		*((size_t*)ptr2) = (size_t)(aligned_ptr - ptr);

		return aligned_ptr;
	}
};

template <unsigned int nByteAlign>
struct AlignedFree{
	void operator()(void *ptr){
		size_t *ptr2=(size_t*)ptr - 1;
		char *p = (char*)ptr;
		p -= (*ptr2);
		free(p);
	}
};

#endif // _MEMORY_ALIGNMENT_H_
