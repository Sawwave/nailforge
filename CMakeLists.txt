cmake_minimum_required(VERSION 3.26)
project(nailforge)

set(CMAKE_CXX_STANDARD 17)


# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -fsanitize=address -Wextra  -Wall") # Debug Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wextra  -Wall -march=native -fno-math-errno -fno-trapping-math")

find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()


add_library(nailforge STATIC        
    src/Alphabet/Alphabet.cpp
    src/PhmmProcessor/PhmmProcessor.cpp
    src/Alphabet/LetterConversion.cpp
    src/SeedExtension/SeedExtension.cpp
    src/StringTree/StringTree.cpp
    src/StringTree/MaxExtensionTable/MaxExtensionTable.cpp
    src/nailforge.cpp)

include_directories(src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/P7HmmReader/build/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/P7HmmReader/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/AvxWindowFmIndex/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/AvxWindowFmIndex/lib/FastaVector/src)

add_subdirectory(AvxWindowFmIndex/)
target_link_libraries(nailforge PRIVATE awfmindex)

add_custom_target(default_target 
    COMMENT "Bundling static libraries..."
    COMMENT "Building p7hmmReader..."
    COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/P7HmmReader -f makefile static    
)

# Custom target for building submodules
add_custom_target(build_awfmi
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/AvxWindowFmIndex/
    COMMENT "Building awfmindex"
)

# set_target_properties(build_awfmi PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_dependencies(nailforge build_awfmi)

