/* SimpleVector.h
 *
 * ABOUT: A primitive version of std::vector which (unlike std::vector) is allowed
 *        to store objects that don't have a constructor, like C arrays (e.g. float[3])
 * NOTES:
 *  - The original purpose of SimpleVector was to design a C++ wrapper around
 *    malloc'd arrays of C arrays (since e.g. std::vector<float[3]> is illegal)
 *  - A SimpleVector ... 
 *    - manages its own memory and stores its size
 *    - does not support move semantics
 *    - has an operator[] for element access
 */

#ifndef SIMPLE_VECTOR_H
#define SIMPLE_VECTOR_H

#include <cstdlib>
#include <exception>
#include <stdexcept>

template<typename T>
class SimpleVector
{
 public:
	
	SimpleVector();
	SimpleVector(const size_t size);
	~SimpleVector();

	// Memory management
	void   resize(const size_t size);
	void   clear();
	size_t size() const { return size_; }

	// Data access
	T&       operator[] (const size_t i)       { return data_[i]; }
	const T& operator[] (const size_t i) const { return data_[i]; }

	T*       data()       { return data_; }
	const T* data() const { return data_; }

 private:
	// Default values guarantee that a new SimpleVector will
	// always be empty and have its data point to null
	T*     data_ = nullptr;
	size_t size_ = 0;
};


//----- Constructors and Destructor -----//

template<typename T>
SimpleVector<T>::SimpleVector() {}

template<typename T>
SimpleVector<T>::SimpleVector(const size_t size)
{
	this->resize(size);
}

template<typename T>
SimpleVector<T>::~SimpleVector() 
{
	if ( data_ != nullptr ) {
		free(data_);
		data_ = nullptr;
		size_ = 0;
	}
}


//----- Memory Management -----//

template<typename T>
void SimpleVector<T>::resize(const size_t new_size) 
{
	// Avoid the ambiguous behavior of realloc when it tries to change the
	// size of the allocated memory block to 0
	if ( new_size == 0 ) {
		this->clear();
		return;
	}

	// Note: when provided a null input ptr, realloc acts just like malloc
	T* temp_data_ptr = (T*) realloc(data_, new_size*sizeof(T));
	size_ = new_size;

	if ( temp_data_ptr != nullptr ) {
		data_         = temp_data_ptr;
		temp_data_ptr = nullptr;
	}
	else {
		// Realloc returned nullptr, so it must have failed
		const char* message = "realloc failed";
		throw std::runtime_error(message);
	}
}

template<typename T>
void SimpleVector<T>::clear()
{
	if ( size_ > 0 ) {
		free(data_);
		data_ = nullptr;
		size_ = 0;
	}
}


//----- Data Access -----//


#endif /* define SIMPLE_VECTOR_H */
