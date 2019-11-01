#include "System.h"

// Implementations of System::get_cache_line_size(), adapted from code
// written by Nick Strupat (https://github.com/NickStrupat/CacheLineSize)

// Linux
#if defined(__linux__)
#include <unistd.h>
std::size_t System::get_cache_line_size()
{
	return static_cast<std::size_t>( sysconf(_SC_LEVEL1_DCACHE_LINESIZE) );
}

// Mac
#elif defined(__APPLE__)
#include <sys/sysctl.h>
std::size_t System::get_cache_line_size()
{
	std::size_t size = 0;
	std::size_t sizeof_size = sizeof(size);
	sysctlbyname("hw.cachelinesize", &size, &sizeof_size, 0, 0);
	return size;
}

// Windows
#elif defined(_WIN32)
#include <stdlib.h>
#include <windows.h>
std::size_t System::get_cache_line_size()
{
	std::size_t size = 0;
	DWORD buffer_size = 0;
	DWORD i = 0;
	SYSTEM_LOGICAL_PROCESSOR_INFORMATION * buffer = 0;

	GetLogicalProcessorInformation(0, &buffer_size);
	buffer = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION *) malloc(buffer_size);
	GetLogicalProcessorInformation(&buffer[0], &buffer_size);

	for (i = 0; i != buffer_size / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION); ++i) {
		if (buffer[i].Relationship == RelationCache && buffer[i].Cache.Level == 1) {
			size = buffer[i].Cache.LineSize;
			break;
		}
	}

	free(buffer);
	return size;
}

// Other
#else
std::size_t System::get_cache_line_size() {
	return 64;  // 64 bytes is a widely-used size, and therefore a reasonable fallback
}
//#error Unable to implement System::get_cache_line_size - Unrecognized platform

#endif // implementations of System::get_cache_line_size()
