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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcompare-debug-second -Wall -Wno-unknown-pragmas")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  add_definitions(-DGENERIC_USE_BOOST_CCMATH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-return-std-move -Wno-nan-infinity-disabled -Wno-missing-template-arg-list-after-template-kw")
endif()

message("-- Build Type: " ${CMAKE_BUILD_TYPE})
message("-- CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

set(EXE_SUFFIX ".exe")
set(DYLIB_SUFFIX "so")
if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
  add_definitions(-DBOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED)
  include_directories(/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include)
  find_program(HOMEBREW_EXECUTABLE brew)
  mark_as_advanced(FORCE HOMEBREW_EXECUTABLE)
  if (HOMEBREW_EXECUTABLE)
    execute_process(COMMAND ${HOMEBREW_EXECUTABLE} --prefix
      OUTPUT_VARIABLE HOMEBREW_INSTALL_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Detected Homebrew with install prefix: ${HOMEBREW_INSTALL_PREFIX}, adding to CMake search paths.")
    list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_INSTALL_PREFIX}")
    set(DYLIB_SUFFIX "dylib")
  endif()
endif()

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
	set(EXE_SUFFIX "")
	set(DYLIB_SUFFIX "dll")
endif()

# Link Boost
find_package(Boost CONFIG REQUIRED COMPONENTS serialization)
include_directories(${Boost_INCLUDE_DIRS})
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

if(DEFINED ENV{EIGEN_PATH})
  include_directories($ENV{EIGEN_PATH})
  message(STATUS "Eigen library path: $ENV{EIGEN_PATH}")
else()
  find_package(Eigen3 REQUIRED)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

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

include_directories(${CMAKE_BINARY_DIR}/include ${PNG_INCLUDE_DIR})

get_directory_property(INC_DIRS INCLUDE_DIRECTORIES)
if(INC_DIRS)
  message(STATUS "Include directories: ${INC_DIRS}")
else()
  message(STATUS "Include directories: <none set>")
endif()

# Build `test` if test.cpp exists in this folder
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test.cpp)
	set(BINARY_NAME test${EXE_SUFFIX})
	add_executable(${BINARY_NAME} test.cpp)
	target_link_libraries(${BINARY_NAME} PRIVATE Boost::serialization Threads::Threads ${PNG_LIBRARY})
endif()

set(BINARY_NAME unit_test${EXE_SUFFIX})
add_executable(${BINARY_NAME} UnitTest.cpp)
target_link_libraries(${BINARY_NAME} PRIVATE Boost::serialization Threads::Threads ${PNG_LIBRARY})