### Test subdirectory CMake

# Use includes copied by top-level CMake
include_directories(${CMAKE_BINARY_DIR}/include ${PNG_INCLUDE_DIR} ${BOOST_INCLUDE_PATH})

# Build `test` if test.cpp exists in this folder
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test.cpp)
	add_executable(test test.cpp)
	target_link_libraries(test PRIVATE Threads::Threads ${PNG_LIBRARY})
	if(TARGET boost_serialization)
		target_link_libraries(test PRIVATE boost_serialization)
	endif()
endif()

add_executable(unit_test UnitTest.cpp)
target_link_libraries(unit_test PRIVATE Threads::Threads ${PNG_LIBRARY})
if(TARGET boost_serialization)
	target_link_libraries(unit_test PRIVATE boost_serialization)
endif()
target_link_libraries(unit_test PRIVATE Threads::Threads ${PNG_LIBRARY})