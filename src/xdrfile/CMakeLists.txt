add_library(indus_xdrfile STATIC
	xdrfile.c
	xdrfile_xtc.c
	xdrfile_trr.c
)

set_target_properties(indus_xdrfile PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(indus_xdrfile PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")

install(TARGETS indus_xdrfile DESTINATION lib)
