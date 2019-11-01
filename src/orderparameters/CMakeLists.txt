### LIB_INDUS ###

# Builds main library: ${LIB_INDUS}
# - TODO: Build both SHARED and OBJECT libraries

# Locate source files
file(GLOB_RECURSE LIB_INDUS_SRC "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

# Declare libraries
# - Don't compile as a static library
#   - It would be preferable to able to create a standalone static library (e.g. '.a') that is then
#     directly linked with various drivers (e.g. normal, testing, ...)
#   - However, trying to do this with GenericFactory leads to awful runtime errors, since
#     static symbols from a static library are frequently optimized out/not initialized in time
# - To mimic a static library, instead compile as an "object" library to allow reuse of '.o' files 
#   when compiling for different drivers
add_library(${LIB_INDUS} OBJECT ${LIB_INDUS_SRC})


function(targetLinkLibraryPreferStatic TARGET SCOPE LIBNAME)
	if(LIB_${LIBNAME}_STATIC)
		target_link_libraries(${TARGET} ${SCOPE} ${LIB_${LIBNAME}_STATIC})
	else()
		target_link_libraries(${TARGET} ${SCOPE} ${LIB_${LIBNAME}})
		target_link_directories(${TARGET} ${SCOPE} "${${LIBNAME}_DIR}/lib")
	endif()
endfunction()


# xdrfile
target_include_directories(${LIB_INDUS} PUBLIC "${XDRFILE_DIR}/include")
# - Prefer static linkage
if(LIB_XDRFILE_STATIC)
	target_link_libraries(${LIB_INDUS} PUBLIC ${LIB_XDRFILE_STATIC})
else()
	# Fall back on dynamic linkage of XDRFILE
	target_link_libraries(${LIB_INDUS} PUBLIC ${LIB_XDRFILE})
	target_link_directories(${LIB_INDUS} PUBLIC "${XDRFILE_DIR}/lib")
endif()

# OpenMP
if(OPENMP_ENABLED)
	target_link_libraries(${LIB_INDUS} PUBLIC OpenMP::OpenMP_CXX)
endif()

# GPTL
if(GPTL_ENABLED)
	target_include_directories(${LIB_INDUS} PUBLIC "${GPTL_DIR}/include")

	# Prefer static linkage
	if(LIB_GPTL_STATIC)
		target_link_libraries(${LIB_INDUS} PUBLIC ${LIB_GPTL_STATIC})
	else()
		# Fall back on dynamic linkage of GPTL
		target_link_libraries(${LIB_INDUS} PUBLIC ${LIB_GPTL})
		target_link_directories(${LIB_INDUS} PUBLIC "${GPTL_DIR}/lib")
	endif()
endif()

# Where the libraries should be installed
install(TARGETS ${LIB_INDUS} DESTINATION lib)
