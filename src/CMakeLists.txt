### LIB_INDUS ###

# This directory contains the core src files
add_subdirectory(orderparameters)


### Driver ###

# Add the executable
add_executable(${INDUS_DRIVER} driver/main.cpp driver/main.h)
target_link_libraries(${INDUS_DRIVER} PUBLIC ${LIB_INDUS})
target_include_directories(${INDUS_DRIVER} PUBLIC ${CMAKE_SOURCE_DIR}/src)

if(APPLE)
	# Run 'codesign' on executables built on Mac
	# - If MPI is enabled, this will prevent the endless firewall prompts about allowing
	#   incoming connections (it will only appear when the program is first run)
	add_custom_command(TARGET ${INDUS_DRIVER}
	                   POST_BUILD
	                   WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
	                   COMMAND codesign --force --deep --sign - ${INDUS_DRIVER}
	                   )
endif()

# Install rules
install(TARGETS ${INDUS_DRIVER} DESTINATION bin)
install(TARGETS ${LIB_INDUS}    DESTINATION lib)
