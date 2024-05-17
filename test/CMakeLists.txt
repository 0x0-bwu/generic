cmake_minimum_required(VERSION 3.14)
project(unit_test.exe)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message("-- Build Type: " ${CMAKE_BUILD_TYPE})

add_compile_options(-Wall)
add_compile_options(-Wno-deprecated-declarations)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-fcompare-debug-second) #hide note message
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	add_definitions(-DGENERIC_USE_BOOST_CCMATH)
  	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-return-std-move")
endif()

set(DYLIB_SUFFIX "so")
# MacOs
# On macOS, add the Homebrew prefix (with appropriate suffixes) to the
# respective HINTS directories (after any user-specified locations).  This
# handles Homebrew installations into non-standard locations (not /usr/local).
# We do not use CMAKE_PREFIX_PATH for this as given the search ordering of
# find_xxx(), doing so would override any user-specified HINTS locations with
# the Homebrew version if it exists.
if (CMAKE_SYSTEM_NAME MATCHES "Darwin")

  	add_definitions(-DBOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED)

  	find_program(HOMEBREW_EXECUTABLE brew)
  	mark_as_advanced(FORCE HOMEBREW_EXECUTABLE)
  	if (HOMEBREW_EXECUTABLE)
    	# Detected a Homebrew install, query for its install prefix.
    	execute_process(COMMAND ${HOMEBREW_EXECUTABLE} --prefix
      	OUTPUT_VARIABLE HOMEBREW_INSTALL_PREFIX
      	OUTPUT_STRIP_TRAILING_WHITESPACE)
    	message(STATUS "Detected Homebrew with install prefix: "
      		"${HOMEBREW_INSTALL_PREFIX}, adding to CMake search paths.")
    	list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_INSTALL_PREFIX}")
    	#set(CMAKE_C_COMPILER "${HOMEBREW_INSTALL_PREFIX}/Cellar/gcc/13.2.0/bin/gcc-13")
    	#set(CMAKE_CXX_COMPILER "${HOMEBREW_INSTALL_PREFIX}/Cellar/gcc/13.2.0/bin/g++-13")
    	set(DYLIB_SUFFIX "dylib")
  	endif()
endif()

# Boost library
set(BOOST_INCLUDE_PATH $ENV{BOOST_PATH}/include)
set(BOOST_LIBRARY_PATH $ENV{BOOST_PATH}/lib)
message(STATUS "BOOST_INCLUDE_PATH: ${BOOST_INCLUDE_PATH}")
message(STATUS "BOOST_LIBRARY_PATH: ${BOOST_LIBRARY_PATH}")

## Ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

## Link Eigen
if(DEFINED ENV{EIGEN_PATH})
	include_directories($ENV{EIGEN_PATH})
else()
	find_package(Eigen3 REQUIRED)
	include_directories(${EIGEN3_INCLUDE_DIR})
endif()

## Link boost_serialization
set(BOOST_SERIALIZATION_LIB_PATH ${BOOST_LIBRARY_PATH}/libboost_serialization.${DYLIB_SUFFIX})
if(NOT EXISTS ${BOOST_SERIALIZATION_LIB_PATH})
	set(BOOST_SERIALIZATION_LIB_PATH ${BOOST_LIBRARY_PATH}/libboost_serialization-mt-x64.${DYLIB_SUFFIX})
endif()
if(EXISTS ${BOOST_SERIALIZATION_LIB_PATH})
	add_library(boost_serialization SHARED IMPORTED)
	add_compile_definitions(GENERIC_BOOST_SERIALIZATION_SUPPORT)
	set_target_properties(boost_serialization PROPERTIES IMPORTED_LOCATION ${BOOST_SERIALIZATION_LIB_PATH})
endif()

find_package(PNG)
if(PNG_FOUND)
add_compile_definitions(GENERIC_BOOST_GIL_IO_PNG_SUPPORT)
endif()

find_package(Threads REQUIRED)

get_filename_component(GENERIC_LIBRARY_DIR ../../ ABSOLUTE)
include_directories(${GENERIC_LIBRARY_DIR} ${PNG_INCLUDE_DIR} ${BOOST_INCLUDE_PATH})

if(EXISTS ${CMAKE_SOURCE_DIR}/test.cpp)
	find_package(PNG REQUIRED)
    add_executable(test.exe test.cpp)
    target_link_libraries(test.exe PRIVATE ${PNG_LIBRARY} Threads::Threads)
	if(EXISTS ${BOOST_SERIALIZATION_LIB_PATH})
    	target_link_libraries(test.exe PRIVATE boost_serialization)
	endif()
endif()
#add_executable(${PROJECT_NAME} UnitTest.cpp)
#target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads ${PNG_LIBRARY})