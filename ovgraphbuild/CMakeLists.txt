cmake_minimum_required (VERSION 2.8)
project (ovgraphbuild)

include_directories ("lib/seqan/include")
set (PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB SOURCE_FILES "${PROJECT_SOURCE_DIR}/*.cpp" "${PROJECT_SOURCE_DIR}/*.h")

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros")

add_executable (ovgraphbuild ${SOURCE_FILES})

find_library (LIBRT_LIB NAMES rt librt)
find_library (PTHREAD_LIB NAMES pthread)

target_link_libraries (ovgraphbuild ${LIBRT_LIB} ${PTHREAD_LIB})

find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories (${ZLIB_INCLUDE_DIRS})
    target_link_libraries (ovgraphbuild ${ZLIB_LIBRARIES})
endif ()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set_property(TARGET ovgraphbuild PROPERTY CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)
