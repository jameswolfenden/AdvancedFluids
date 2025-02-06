cmake_minimum_required(VERSION 3.5.0)
project(fluids VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE RelWithDebInfo)

option(ENABLE_CLANG_TIDY "Enable clang-tidy checks during build" OFF)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

include_directories(${INCLUDE_DIR})

# Collect all source files except main.cpp
file(GLOB_RECURSE SOURCES 
    "${SRC_DIR}/*.cpp"
)
list(FILTER SOURCES EXCLUDE REGEX ".*/main\\.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*/main_old\\.cpp$")

# Print found sources
message(STATUS "Found source files:")
foreach(SOURCE ${SOURCES})
    message(STATUS "  ${SOURCE}")
endforeach()

if(ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES clang-tidy)
    if(CLANG_TIDY_EXE)
        message(STATUS "Found clang-tidy: ${CLANG_TIDY_EXE}")
        set(CMAKE_CXX_CLANG_TIDY 
            ${CLANG_TIDY_EXE};
            -header-filter=^${CMAKE_SOURCE_DIR}/.*;
            -config-file=${CMAKE_SOURCE_DIR}/.clang-tidy
        )
    else()
        message(FATAL_ERROR "clang-tidy requested but not found")
    endif()
endif()

# Create library from sources
add_library(fluids_lib STATIC ${SOURCES})

# Create main executable
add_executable(fluids ${SRC_DIR}/main.cpp)

# Find and link HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)
if (HDF5_FOUND)
    target_include_directories(fluids_lib PUBLIC ${HDF5_INCLUDE_DIRS} ${HDF5_CXX_INCLUDE_DIRS})
    target_link_libraries(fluids_lib PUBLIC ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})
endif()

# Find and link OpenMP
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    target_link_libraries(fluids_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Include google test
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
)
FetchContent_MakeAvailable(googletest)

# Link executable to library (dependencies will propagate)
target_link_libraries(fluids PRIVATE fluids_lib)

# Add compiler warnings
if(MSVC)
    target_compile_options(fluids_lib PRIVATE /W4)
    target_compile_options(fluids PRIVATE /W4)
    target_compile_options(fluids_lib PRIVATE /O2 /fp:fast /arch:AVX2)
else()
    target_compile_options(fluids_lib PRIVATE -Wall -Wextra)
    target_compile_options(fluids PRIVATE -Wall -Wextra)
    target_compile_options(fluids_lib PRIVATE 
    -O3                  # Aggressive optimisation
    -march=native        # Optimize for current CPU
    -ftree-vectorize     # Enable vectorisation (usually on with -O3)
    -ffast-math          # Allow math optimisations
    -fopt-info-vec       # Report vectorisation info (GCC)
    )
endif()

# Add tests if they exist
if(EXISTS "${PROJECT_SOURCE_DIR}/tests/CMakeLists.txt")
    enable_testing()
    add_subdirectory(tests)
endif()

# Add benchmarks if they exist
if(EXISTS "${PROJECT_SOURCE_DIR}/benchmarks/CMakeLists.txt")
    add_subdirectory(benchmarks)
endif()

# Print status information
message(STATUS "HDF5 libs: ${HDF5_LIBRARIES}")
message(STATUS "HDF5 includes: ${HDF5_INCLUDE_DIRS}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
