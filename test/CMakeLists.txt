cmake_minimum_required(VERSION 3.26)
project(nailforgeBenchmark)

set(CMAKE_CXX_STANDARD 20)

add_executable(nailforgeBenchmark
    benchmark.cpp
)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O0 -g -fsanitize=address  -Wextra -march=native")

target_include_directories(nailforgeBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src/)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../P7HmmReader/build/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../P7HmmReader/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../AvxWindowFmIndex/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../AvxWindowFmIndex/lib/FastaVector/src)

add_library(nailforge STATIC IMPORTED)
set_target_properties(nailforge PROPERTIES IMPORTED_LOCATION ../libnailforge.a)

add_library(AwFmIndexLib STATIC IMPORTED)
set_target_properties(AwFmIndexLib PROPERTIES IMPORTED_LOCATION ../AvxWindowFmIndex/build/libawfmindex_static.a)

add_library(P7HmmLib STATIC IMPORTED)
set_target_properties(P7HmmLib PROPERTIES IMPORTED_LOCATION ../P7HmmReader/build/lib/libP7HmmReader.a)

add_library(FastaVectorLib STATIC IMPORTED)
set_target_properties(FastaVectorLib PROPERTIES IMPORTED_LOCATION ../AvxWindowFmIndex/lib/FastaVector/build/libfastavector_static.a)

add_library(libdivsufsort SHARED IMPORTED)
set_target_properties(libdivsufsort PROPERTIES IMPORTED_LOCATION ../AvxWindowFmIndex/lib/libdivsufsort/lib/build/libdivsufsort64.a)
target_link_libraries(nailforgeBenchmark PRIVATE nailforge )
target_link_libraries(nailforgeBenchmark PRIVATE AwFmIndexLib )
target_link_libraries(nailforgeBenchmark PRIVATE libdivsufsort )
target_link_libraries(nailforgeBenchmark PRIVATE FastaVectorLib  )
target_link_libraries(nailforgeBenchmark PRIVATE P7HmmLib)

include_directories(src)


find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

