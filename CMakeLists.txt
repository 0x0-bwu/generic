cmake_minimum_required(VERSION 3.14...4.24)

project(generic)

# Remove previously-copied project include (if any), then copy headers
if(EXISTS "${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}")
  file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}")
  message(STATUS "Removed existing ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}")
endif()

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/include")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/src" DESTINATION "${CMAKE_BINARY_DIR}/include")
file(RENAME "${CMAKE_BINARY_DIR}/include/src" "${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}")

# Option to control building tests
option(BUILD_GENERIC_TEST "Build unit test" OFF)
if(BUILD_GENERIC_TEST)
  message(STATUS "BUILD_GENERIC_TEST=ON: adding test subdirectory")
  add_subdirectory(test)
else()
  message(STATUS "BUILD_GENERIC_TEST=OFF: skipping test subdirectory (set -DBUILD_GENERIC_TEST=ON to enable)")
endif()
