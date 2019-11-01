// System - functions and static wrappers for querying system info

#ifndef SYSTEM_H
#define SYSTEM_H

#include <cstddef>  // defines std::size_t

class System {
 public:
	// Returns the cache line size (in bytes)
	static std::size_t get_cache_line_size();
};

#endif // SYSTEM_H
