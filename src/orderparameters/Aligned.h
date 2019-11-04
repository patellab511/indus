// Aligned.h - aligned dynamic memory allocation
// - Provides general, C-style wrappers around system-specified functions that 
//   allocate/deallocate aligned memory
// - Also provides a std::allocator-style object that can be used with other
//   C++ objects in the STL, such as std::vector
//
// - The code here is based on the folloing prior works
//   - 'Mallocator', by Stephen T. Lavavej
//      - Original:
//      - Updated: slides presented at CppCon 2014 (under https://github.com/CppCon/CppCon2014/blob/master/Presentations)
// 
//   - 'Aligned', by Nick Strupat
//     - GitHub: https://github.com/NickStrupat/Aligned
//
// - Misc. notes
//   - 'construct()' and 'destruct()' are not implemented, as they are redundant in C++11
//
#ifndef ALIGNED_H
#define ALIGNED_H

// The following headers are required for all allocators.
#include <cstddef>    // size_t and ptrdiff_t
#include <exception>
#include <memory>     // allocator_traits
#include <new>        // placement new, bad_alloc, bad_array_new_length
#include <stdexcept>  // std::length_error

#include <cstdlib>    // malloc() and free()
#include <iostream>   // std::cout
#include <limits>
#include <ostream>    // std::endl

#include "System.h"

namespace Aligned {

using size_type = std::size_t;

// Aligned malloc
void* malloc(const size_type size, const size_type alignment)
{
	// Require alignment/sizeof(void*) = 2^n for some positive integer 'n'
	size_type n = alignment/(sizeof(void*));
	bool is_power_of_2 = ( (n > 0) and ((n & (n-1)) == 0) );
	if ( not is_power_of_2 ) {
		throw std::runtime_error("Error in Aligned::malloc - alignment/sizeof(void*) must be a power of 2");
	}

	if ( n == 0 ) {
		return nullptr;
	}

	// Allocate memory
	int ret;
	void* mem_ptr = nullptr;
#if defined(__GNUC__) || defined(__APPLE__)
	ret = posix_memalign(&mem_ptr, alignment, size);
	if ( ret != 0 ) {
		if ( mem_ptr != nullptr ) { free(mem_ptr); }
		throw std::runtime_error("posix_memalign: ret is nonzero");
	}
#elif defined(__intel__)
	mem_ptr = _mm_alloc(size, alignment);
#elif defined(_WIN32)
	mem_ptr = _aligned_malloc(size, alignment);
#else
	// TODO As a fallback: allocate memory_needed+alignment bytes,
	//  and find an address that matches the criteria
	static_assert(false, "Unsupported aligned malloc");
#endif

	return mem_ptr;
}


// Aligned free
void free(void* mem_ptr) noexcept
{
#if defined(__GNUC__) || defined(__APPLE__)
	// Memory allocated by 'posix_memalign()' is deallocated by the usual 'free()'
	::free(mem_ptr);
#elif defined(__intel__)
	_mm_free(mem_ptr);
#elif defined(_WIN32)
	_aligned_free(mem_ptr);
#else
	::free(mem_ptr);
#endif
}


// STL-compatible allocator for aligned memory
// - If template parameter 'Alignment' is zero or unspecified,
//   data is aligned according to the size of the cache line
//   (which must be queried at runtime)
template <typename T, size_type Alignment = 0>
//class Allocator {
class Allocator {
 public:
	// The following will be the same for virtually all allocators.
	using value_type      = T;
	using pointer         = T*;
	using const_pointer   = const T*;
	using reference       = T&;
	using const_reference = const T&;
	using difference_type = std::ptrdiff_t;

	// When rebinding, use the same alignment
	template <typename U> 
	struct rebind {
		using other = Allocator<U>;
	};

	// Returns true if and only if storage allocated from *this
	// can be deallocated from other, and vice versa.
	bool operator==(const Allocator& other) const noexcept {
		return true;
	}

	bool operator!=(const Allocator& other) const noexcept {
		// Given the implementation of operator== above, this always returns false
		return (not ( *this == other ));  
	}

	// Default constructor: use template argument as alignment
	Allocator() {
		set_alignment(Alignment);
	}

	// Copy constructor: copy alignment from other allocator
	Allocator(const Allocator& other) {
		set_alignment(other.alignment_);
	}

	// Rebinding constructor, for making allocators of different types
	template <typename U> 
	Allocator(const Allocator<U>& other) {
		set_alignment(other.alignment_);
	}

	~Allocator() { }

	pointer allocate(const size_type n) const {
		void* const pv = Aligned::malloc(n * sizeof(value_type), alignment_);

		if ( pv == nullptr ) {
			// Allocation failed
			throw std::bad_alloc();
		}

		return static_cast<pointer>(pv);
	}

	void deallocate(pointer const p, const size_type n) const noexcept {
		Aligned::free( static_cast<void*>(p) );
	}

	size_type get_alignment() const noexcept {
		return alignment_;
	}

 private:
	size_type alignment_ = 0;

	void set_alignment(const size_type alignment) {
		if ( alignment != 0 ) {
			alignment_ = alignment;
		}
		else {
			alignment_ = System::get_cache_line_size();
		}
	}

	// Allocators are not required to be assignable, so all allocators should have a private 
	// unimplemented assignment operator
	Allocator& operator=(const Allocator&);
};

} // end namespace Aligned
#endif // ifndef ALIGNED_H
