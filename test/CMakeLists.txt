# CMake setup
cmake_minimum_required(VERSION 3.14)

# Project name
project(generic)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(NOT CMAKE_BUILD_TYPE)
  	set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

# Compile options
option(ENABLE_ASSERTION "Enable assertion" OFF)
option(ENABLE_EXCEPTION "Enable exception" ON)
if(NOT ENABLE_ASSERTION)
	add_compile_definitions(GENERIC_NO_ASSERTION)	
endif()
if(NOT ENABLE_EXCEPTION)
	add_compile_definitions(GENERIC_NO_EXCEPTION)
endif()
message(STATUS "Enable assertion: ${ENABLE_ASSERTION}")
message(STATUS "Enable exception: ${ENABLE_EXCEPTION}")

option(ENABLE_PROFILING "Enable profiling" OFF)
if (ENABLE_PROFILING)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
	message("-- Enable GNU Profiling")
endif()

option(ENABLE_ASAN "Enable address sanitizer" OFF)
option(ENABLE_TSAN "Enable thread sanitizer" OFF)

if (ENABLE_ASAN AND ENABLE_TSAN) 
	message(FATAL_ERROR "Only one of ENABLE_ASAN, ENABLE_TSAN can be enabled at the same time.")
endif()

if (ENABLE_ASAN)
	add_compile_options(-fsanitize=address)
	add_link_options(-fsanitize=address)
endif()

if (ENABLE_TSAN)
	add_compile_options(-fsanitize=thread)
	add_link_options(-fsanitize=thread)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -ffast-math -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -ffast-math -Wno-deprecated-declarations")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcompare-debug-second -Wall -Wno-unknown-pragmas") #hide note message and clang pragma warning
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	add_definitions(-DGENERIC_USE_BOOST_CCMATH)
  	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-return-std-move -Wno-nan-infinity-disabled -Wno-missing-template-arg-list-after-template-kw")
endif()

message("-- Build Type: " ${CMAKE_BUILD_TYPE})
message("-- CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

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
	include_directories(/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include)
	
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

# Windows / MSYS2 support: prefer .dll and use common MSYS2 MinGW prefixes
if (CMAKE_SYSTEM_NAME MATCHES "Windows" OR WIN32)
	set(DYLIB_SUFFIX "dll")
	if(DEFINED ENV{MSYSTEM})
		message(STATUS "Detected MSYS2 environment: $ENV{MSYSTEM}")
		if("$ENV{MSYSTEM}" STREQUAL "MINGW64")
			set(MSYS2_MINGW_PREFIX "/mingw64")
		elseif("$ENV{MSYSTEM}" STREQUAL "MINGW32")
			set(MSYS2_MINGW_PREFIX "/mingw32")
		else()
			set(MSYS2_MINGW_PREFIX "/mingw64")
		endif()
		if(BOOST_INCLUDE_PATH STREQUAL "" AND BOOST_LIBRARY_PATH STREQUAL "")
			set(BOOST_INCLUDE_PATH ${MSYS2_MINGW_PREFIX}/include)
			set(BOOST_LIBRARY_PATH ${MSYS2_MINGW_PREFIX}/lib)
			message(STATUS "Using MSYS2 Boost paths: ${BOOST_INCLUDE_PATH}, ${BOOST_LIBRARY_PATH}")
		endif()
  	endif()
endif()

# Boost library
if(DEFINED ENV{BOOST_PATH})
set(BOOST_INCLUDE_PATH $ENV{BOOST_PATH}/include)
set(BOOST_LIBRARY_PATH $ENV{BOOST_PATH}/lib)
else()
	set(BOOST_INCLUDE_PATH "")
	set(BOOST_LIBRARY_PATH "")
endif()
message(STATUS "Boost include path: ${BOOST_INCLUDE_PATH}")
message(STATUS "Boost library path: ${BOOST_LIBRARY_PATH}")

## Ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

## Link boost_serialization
set(BOOST_SERIALIZATION_LIB_PATH "")
set(BOOST_SERIALIZATION_CANDS
	${BOOST_LIBRARY_PATH}/libboost_serialization.${DYLIB_SUFFIX}
	${BOOST_LIBRARY_PATH}/libboost_serialization-mt.${DYLIB_SUFFIX}
	${BOOST_LIBRARY_PATH}/libboost_serialization-mt-x64.${DYLIB_SUFFIX}
	${BOOST_LIBRARY_PATH}/boost_serialization.${DYLIB_SUFFIX}
)
foreach(CANDIDATE ${BOOST_SERIALIZATION_CANDS})
	if(EXISTS ${CANDIDATE})
		set(BOOST_SERIALIZATION_LIB_PATH ${CANDIDATE})
		break()
endif()
endforeach()
if(BOOST_SERIALIZATION_LIB_PATH)
	add_library(boost_serialization SHARED IMPORTED)
	add_compile_definitions(GENERIC_BOOST_SERIALIZATION_SUPPORT)
	set_target_properties(boost_serialization PROPERTIES IMPORTED_LOCATION ${BOOST_SERIALIZATION_LIB_PATH})
endif()

# Link Eigen
if(DEFINED ENV{EIGEN_PATH})
	include_directories($ENV{EIGEN_PATH})
	message(STATUS "Eigen library path: $ENV{EIGEN_PATH}")
else()
	find_package(Eigen3 REQUIRED)
	include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# Find package
find_package(PNG)
if(PNG_FOUND)
	add_compile_definitions(GENERIC_BOOST_GIL_IO_PNG_SUPPORT)
endif()

find_package(Threads REQUIRED)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "FOUND OMP: ${OpenMP_CXX_LIBRARY}")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Copy project headers to include under build dir 
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/include")
file(COPY "${CMAKE_SOURCE_DIR}/../src" DESTINATION "${CMAKE_BINARY_DIR}/include")
file(RENAME "${CMAKE_BINARY_DIR}/include/src" "${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}")

include_directories(${CMAKE_BINARY_DIR}/include ${PNG_INCLUDE_DIR} ${BOOST_INCLUDE_PATH})

# Print configured include directories so it's easy to debug include search paths
get_directory_property(INC_DIRS INCLUDE_DIRECTORIES)
if(INC_DIRS)
	message(STATUS "Include directories: ${INC_DIRS}")
else()
	message(STATUS "Include directories: <none set>")
endif()

# Add Executable
if(EXISTS ${CMAKE_SOURCE_DIR}/test.cpp)
    add_executable(test test.cpp)
    target_link_libraries(test PRIVATE ${PNG_LIBRARY} Threads::Threads)
	if(EXISTS ${BOOST_SERIALIZATION_LIB_PATH})
    	target_link_libraries(test PRIVATE boost_serialization)
	endif()
endif()

add_executable(unit_test UnitTest.cpp)
target_link_libraries(unit_test PRIVATE Threads::Threads ${PNG_LIBRARY})